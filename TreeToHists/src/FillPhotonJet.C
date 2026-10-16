#include "ScientificProducerEvidence.h"
#include "SelectionViewContract.h"
#include "PurityInputs.h"
#include "CandidateCache.h"
#include "MeasurementAccumulator.h"
#include "EmsComposition.h"
#include "LeadingAssociation.h"
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
#include <set>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/* Fill the photon + recoil-jet histograms from eventTree v1 files. */

/* root -l -q 'TreeToHists/src/FillPhotonJet.C("a.eventtree.root b.eventtree.root", "out.root", "nominal")' */

/* Analysis events: evt_record == 0 && evt_nodes_ready == 1 && vtx_valid == 1 && an accepted producer class
   (acceptedEventClass in PhotonJet.h: evt_terminal_status 0, or 4 for a sample stitched here by truth-jet pT).
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

/* A cut-only view reuses this source's reader, scores and subtracted witnesses.
   Its sealed profile must be validated against the parent's invariant body by
   NativeSelectionViews.h before entering the reader. */
struct SelectionMeasurementView {
  Profile profile;
  Histograms* histograms=nullptr;
  EventTable* records=nullptr;
  Histograms* timing=nullptr;
  EventTable* timingRecords=nullptr;
  purity_inputs::Book* purity=nullptr;
  double expectedStoredGap=std::numeric_limits<double>::quiet_NaN();
  Long64_t originalEntries=0,scienceReaders=0;
};

/* Two exact pp DATA JES products share original arrays and accepted regions.
   Factors are absolute; no selected-record replay or second calibration. */
struct DataJesView {
  double scale;
  Long64_t originalEntries=0,scienceReaders=0;
  Histograms histograms, timing;
  EventTable records, timingRecords;
  explicit DataJesView(double value,const std::vector<CentralityClass>& classes)
      : scale(value),histograms(classes),timing(classes) {}
};


inline void fillFile(const std::string& path, int fileIndex, const Profile& requested,
                     const std::vector<CentralityClass>& classes, Histograms& h, EventTable& table,
                     const CentralityWeights& weights = {}, const std::string& dataTrigger = "",
                     Long64_t progressOffset = 0, candidate_cache::Writer* cacheWriter = nullptr, Histograms* timeHists = nullptr,
                     EventTable* timeTable = nullptr,const response::Variation* variation=nullptr,
                     const char* emsContractText="",const std::function<std::string(const std::string&)>& digest={},
                     const std::vector<DataJesView*>& dataJes={},
                     const std::vector<SelectionMeasurementView*>& selectionViews={},
                     purity_inputs::Book* nominalPurity=nullptr) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
  const Profile profile = validateSelectionInput(*file, requested);
  if((nominalPurity || !selectionViews.empty()) && isSimulation(*file))throw std::runtime_error("cut measurement fanout uses DATA reader; SIM uses response reader");
  if(nominalPurity && nominalPurity->family!="data")throw std::runtime_error("DATA purity book family differs");
  std::vector<Profile> selectionProfiles;
  for(auto* view:selectionViews){
    if(!view || !view->histograms || !view->records)throw std::runtime_error("selection measurement view missing accumulator");
    if(view->purity && view->purity->family!="data")throw std::runtime_error("DATA view purity book family differs");
    selection_views::compatible(profile,view->profile,view->expectedStoredGap);
    selectionProfiles.push_back(validateSelectionInput(*file,view->profile));++view->scienceReaders;
  }
  if((nominalPurity || !selectionViews.empty()) && variation && response::variationJson(*variation)!=response::variationJson(response::Variation{}))throw std::runtime_error("selection views cannot combine energy effects");
  if(!dataJes.empty() && (dataJes.size()!=2 || !dataJes[0] || !dataJes[1] || dataJes[0]->scale!=1./.907 || dataJes[1]->scale!=1./.923))
    throw std::runtime_error("fused DATA JES requires its two exact positive finite absolute scales");
  if(!dataJes.empty() && (profile.system!=System::pp || isSimulation(*file)))
    throw std::runtime_error("fused DATA JES accepts original pp DATA only");
  if (fileSystem(*file) != profile.system) throw std::runtime_error("profile system differs from " + path);
  TTreeReader reader("eventTree", file.get());
  for(auto* view:dataJes)++view->scienceReaders;
  applyEntryRange(reader);
  const VertexEventWeight vertexWeight(*file,reader);
  const PpSimulationGate ppSelection(*file,reader);
  const AuAuTreeWeightCopies treeCopies(*file,reader,weights);
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
  TTreeReaderArray<double> originalEt(reader, "pho_et"), eta(reader, "pho_eta"), phi(reader, "pho_phi");
  PhotonScores originalScore(*file,reader,profile);PhotonPreselection originalPresel(reader);
  ems_native::Overlay<TTreeReaderArray<double>,double> et(originalEt);
  ems_native::Overlay<PhotonScores,double> score(originalScore);
  ems_native::Overlay<PhotonPreselection,int> presel(originalPresel);
  std::unique_ptr<TTreeReaderArray<double>> rawTime;
  if(timeHists)rawTime=std::make_unique<TTreeReaderArray<double>>(reader,"pho_time_ns");
  TTreeReaderArray<int> originalState(reader, profile.stateBranch.c_str());
  TTreeReaderArray<double> originalCone(reader, "isow_cone_sum"), originalThreshold(reader, "isow_threshold"),
      originalSideband(reader, "isow_sideband_threshold");
  ems_native::Overlay<TTreeReaderArray<int>,int> state(originalState);
  ems_native::Overlay<TTreeReaderArray<double>,double> cone(originalCone),threshold(originalThreshold),sideband(originalSideband);
  TTreeReaderArray<int> method(reader, "isow_subtraction_method"), role(reader, "isow_reconstructed_or_truth");
  TTreeReaderArray<int> idClass(reader, "pho_id_class"), isoClass(reader, "pho_iso_class");
  TTreeReaderArray<double> storedJetPt(reader, profile.jetBranch("jet", "pt").c_str()), jetEta(reader, profile.jetBranch("jet", "eta").c_str()),
      jetPhi(reader, profile.jetBranch("jet", "phi").c_str());
  ems_native::Overlay<TTreeReaderArray<double>,double> selectedJetPt(storedJetPt);
  const response::Variation chosen=variation?*variation:response::Variation{};
  std::unique_ptr<ems_native::Runtime> emsRuntime;std::unique_ptr<energy_selection::NativeInputs> emsInputs;
  std::unique_ptr<TTreeReaderArray<double>> emFraction;
  if(ems_native::active(profile,chosen)){
    emsRuntime=std::make_unique<ems_native::Runtime>(profile,ems_native::parseContract(emsContractText),digest);
    emsInputs=std::make_unique<energy_selection::NativeInputs>(*file->Get<TTree>("eventTree"),reader,profile,originalScore);
    emFraction=std::make_unique<TTreeReaderArray<double>>(reader,profile.jetBranch("jet","emcal_frac").c_str());
  }
  const double jetScale = isSimulation(*file) ? 1.0 : configuredDataJetScale;  /* the in-situ scale of DATA jets */
  const OwnershipGate owner(reader);  /* a stitched sample's events outside its window belong to another sample */

  const auto scientificContext=[&](const Profile& p,double scale,bool timing) {
    return photonjet_science::measurementContext(*file,p,weights,
        photonjet_science::Json::parse(response::variationJson(chosen)),photonjet_science::Json::object(),scale,timing,dataTrigger,
        originalScore.scientificFunctions());
  };
  h.scientificContext(scientificContext(profile,jetScale,false));
  table.scientificContext(scientificContext(profile,jetScale,false));
  if(timeTable)timeTable->scientificContext(scientificContext(profile,jetScale,true));
  if(timeHists)timeHists->scientificContext(scientificContext(profile,jetScale,true));
  for(auto* view:dataJes) {
    view->records.scientificContext(scientificContext(profile,view->scale,false));
    view->timingRecords.scientificContext(scientificContext(profile,view->scale,true));
    view->histograms.scientificContext(scientificContext(profile,view->scale,false));
    view->timing.scientificContext(scientificContext(profile,view->scale,true));
  }
  for(std::size_t i=0;i<selectionViews.size();++i) {
    selectionViews[i]->records->scientificContext(scientificContext(selectionProfiles[i],jetScale,false));
    if(selectionViews[i]->timingRecords)selectionViews[i]->timingRecords->scientificContext(scientificContext(selectionProfiles[i],jetScale,true));
    selectionViews[i]->histograms->scientificContext(scientificContext(selectionProfiles[i],jetScale,false));
    if(selectionViews[i]->timing)selectionViews[i]->timing->scientificContext(scientificContext(selectionProfiles[i],jetScale,true));
  }

  std::unique_ptr<candidate_cache::Source> cache;
  if(cacheWriter)cacheWriter->scientificContext(fileIndex,scientificContext(profile,jetScale,false));
  if(cacheWriter)cache=std::make_unique<candidate_cache::Source>(*cacheWriter,*static_cast<TTree*>(file->Get("eventTree")),reader,profile.simulation);
  while (reader.Next()) {
    for(auto* view:dataJes)++view->originalEntries;
    for(auto* view:selectionViews)++view->originalEntries;
    et.clear();score.clear();presel.clear();state.clear();cone.clear();threshold.clear();sideband.clear();selectedJetPt.clear();
    if(cache)cache->observeSource(reader.GetCurrentEntry());
    treeCopies.check();
    reportProgress(progressOffset+reader.GetCurrentEntry()+1);
    validateSelectionEvent(profile, *centrality, originalEt, originalScore, originalState, originalCone, originalThreshold, originalSideband, method, role, idClass, isoClass);
    if (*record != 0 || *ready != 1 || *vtxValid != 1 || !acceptedEventClass(*terminal, *centrality)) continue;
    if (!ppSelection.accepts()) continue;
    if (!owner.owned()) continue;
    if(profile.configured && (!std::isfinite(*vz) || std::fabs(*vz)>=profile.vertexAbsZMax)) continue;
    if(!dataSelection.accepts(*run,scaled?**scaled:0,mb?**mb:0)) continue;
    if(!dataTrigger.empty() && !(**scaled & triggerMask)) continue;
    const candidate_cache::WeightFactors factors{*weightValue,treeCopies.baseWeight(ppSelection.weight(*weightValue)),vertexWeight(),chosen.centralityOff ? 1.0 : weights(*centrality)};
    const double product = factors.base * factors.vertex * factors.centrality;
    if (strictReduction && !std::isfinite(product)) throw std::runtime_error("nonfinite event weight in " + path);
    const double w = std::isfinite(product) ? product : 1.0;
    if(cache)cache->append(*run,*seq,reader.GetCurrentEntry(),fileIndex,*centrality,*vz,w,[&](size_t i){return score[i];},factors,[&](size_t i){
      const auto decision=photonClass(profile,i,et[i],*centrality,score,state,cone,threshold,sideband,method,role,idClass,isoClass);
      auto values=originalScore.cacheInputs(i);
      values.insert({{"applied_id_class",decision.tight?kIdTight:decision.nonTight?kIdNonTight:decision.idCode},
        {"applied_iso_class",decision.isolated?kIsoIsolated:decision.nonIsolated?kIsoNonIsolated:decision.isoCode},
        {"applied_region",decision.region()},{"applied_preselection",presel[i]},
        {"applied_kinematics",photonKinematics(profile,*vz,et[i],eta[i],phi[i])},
        {"analysis_isolation",analysisIsolation(profile,cone[2*i+profile.isolationRadiusIndex])}});
      return values;
    });
    if(emsRuntime){
      const auto selected=emsRuntime->evaluate(*emsInputs,originalEt,*vz,*centrality,chosen.emScale);
      et.assign(selected.et);score.assign(selected.score);presel.assign(selected.preselection);state.assign(selected.state);
      cone.assign(selected.cone);threshold.assign(selected.threshold);sideband.assign(selected.sideband);
      response::Varied jets(storedJetPt,1.,0.,chosen.seed,0,reader.GetCurrentEntry(),1);
      jets.scaleEmFraction(*emFraction,chosen.emScale);selectedJetPt.assign(jets);
    }
    int k = 0;  /* the event's centrality class (classes are disjoint) */
    if (profile.system == System::auau) {
      while (k < static_cast<int>(classes.size()) && !(classes[k].low <= *centrality && *centrality < classes[k].high)) ++k;
      if (k == static_cast<int>(classes.size())) continue;
    }
    const CentralityClass& c = classes[k];
    h.h1("events_" + c.name).Fill(0.5, w);
    if(timeHists)timeHists->h1("events_"+c.name).Fill(.5,w);
    for(auto* view:dataJes){view->histograms.h1("events_"+c.name).Fill(.5,w);view->timing.h1("events_"+c.name).Fill(.5,w);}
    for(auto* view:selectionViews){view->histograms->h1("events_"+c.name).Fill(.5,w);
      if(view->timing)view->timing->h1("events_"+c.name).Fill(.5,w);}
    const size_t nPhoton = et.GetSize();
    if (nPhoton == 0) continue;
    if (cone.GetSize() != 2 * nPhoton) throw std::runtime_error("isolation witnesses are not 2 per photon in " + path);
    if (score.GetSize() != nPhoton || (!profile.optionalID && state.GetSize() != nPhoton))
      throw std::runtime_error("photon-ID branches pho_bdt_score, pho_bdt_state are not filled in " + path);
    if (profile.augmentedClasses && (idClass.GetSize() != nPhoton || isoClass.GetSize() != nPhoton))
      throw std::runtime_error("the augmented classes pho_id_class, pho_iso_class are not filled in " + path);
    auto accumulate = [&](Histograms& target, EventTable* saved, const std::vector<int>& accepted, bool timing, double scale, const std::array<long,4>* leaders=nullptr) {
      std::vector<int> visible(nPhoton,1);
      if(timing)for(size_t i=0;i<nPhoton;++i)visible[i]=dataTimeCut((*rawTime)[i]);
      accumulateMeasurement(profile,c,*run,*seq,reader.GetCurrentEntry(),fileIndex,*centrality,*vz,w,
          target,saved,et,eta,phi,score,accepted,visible,state,cone,threshold,sideband,method,role,
          idClass,isoClass,selectedJetPt,jetEta,jetPhi,scale,leaders);
    };
    std::vector<int> accepted(nPhoton);for(size_t i=0;i<nPhoton;++i)accepted[i]=presel[i];
    std::array<long,4> leaders;
    if(!dataJes.empty())leaders=regionLeaders(profile,*vz,*centrality,nPhoton,et,eta,phi,score,accepted,state,cone,threshold,sideband,method,role,idClass,isoClass);
    const auto* shared=dataJes.empty()?nullptr:&leaders;
    accumulate(h,&table,accepted,false,jetScale,shared);
    if(nominalPurity)nominalPurity->accumulate(profile,c,*vz,*centrality,w,originalEt,eta,phi,originalScore,
      accepted,originalState,originalCone,originalThreshold,originalSideband,method,role,idClass,isoClass,storedJetPt,jetEta,jetPhi,
      std::vector<int>(nPhoton,0),jetScale);
    for(size_t v=0;v<selectionViews.size();++v){
      auto& target=*selectionViews[v];const auto& variedProfile=selectionProfiles[v];
      selection_views::validateStoredGap(profile,target.expectedStoredGap,nPhoton,originalCone,originalThreshold,originalSideband,method,role);
      accumulateMeasurement(variedProfile,c,*run,*seq,reader.GetCurrentEntry(),fileIndex,*centrality,*vz,w,
        *target.histograms,target.records,et,eta,phi,score,accepted,std::vector<int>(nPhoton,1),state,
        cone,threshold,sideband,method,role,idClass,isoClass,selectedJetPt,jetEta,jetPhi,jetScale);
      if(target.purity)target.purity->accumulate(variedProfile,c,*vz,*centrality,w,originalEt,eta,phi,originalScore,
        accepted,originalState,originalCone,originalThreshold,originalSideband,method,role,idClass,isoClass,storedJetPt,jetEta,jetPhi,
        std::vector<int>(nPhoton,0),jetScale);
      if(target.timing){
        if(!rawTime || !target.timingRecords || rawTime->GetSize()!=nPhoton)throw std::runtime_error("selection timing needs original timing reader");
        auto timed=accepted;std::vector<int> visible(nPhoton);
        for(size_t i=0;i<nPhoton;++i){visible[i]=dataTimeCut((*rawTime)[i]);timed[i]=timed[i]&&visible[i];}
        accumulateMeasurement(variedProfile,c,*run,*seq,reader.GetCurrentEntry(),fileIndex,*centrality,*vz,w,
          *target.timing,target.timingRecords,et,eta,phi,score,timed,visible,state,cone,threshold,sideband,method,role,
          idClass,isoClass,selectedJetPt,jetEta,jetPhi,jetScale);
      }
    }
    for(auto* view:dataJes)accumulate(view->histograms,&view->records,accepted,false,view->scale,shared);
    if(timeHists){
      if(rawTime->GetSize()!=nPhoton)throw std::runtime_error("DATA timing array size differs");
      for(size_t i=0;i<nPhoton;++i)accepted[i]=accepted[i]&&dataTimeCut((*rawTime)[i]);
      if(!dataJes.empty())leaders=regionLeaders(profile,*vz,*centrality,nPhoton,et,eta,phi,score,accepted,state,cone,threshold,sideband,method,role,idClass,isoClass);
      accumulate(*timeHists,timeTable,accepted,true,jetScale,shared);
      for(auto* view:dataJes)accumulate(view->timing,&view->timingRecords,accepted,true,view->scale,shared);
    }
  }
  /* a missing branch or a read error ends the loop like the end of the tree; only the end is accepted */
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd && reader.GetEntryStatus() != TTreeReader::kEntryNotFound)
    throw std::runtime_error("eventTree read stopped before the end of " + path);
  if(cache)cache->finishSource(fileIndex);
}

}  /* namespace photonjet */

/* inputs: whitespace-separated eventTree paths. profile: "nominal" (the analysis), a variant of it (PhotonJet.h
   namedProfile), or "historical" (the reference reproduction only). */
/* Region-leading photons and their recoil jets retain event references
   alongside histogram cells, preserving correlations when FinalAnalysis
   draws one bootstrap multiplier per event. */
int FillPhotonJet(const char* inputs, const char* output, const char* profileName = "nominal", const char* selection = "",
                  const char* centralityWeightText = "",const char* variationText="",const char* emsContractText="",
                  const std::function<std::string(const std::string&)>& digest={},
                  const char* dataJesUp=nullptr,const char* dataJesDown=nullptr) {
  using namespace photonjet;
  const auto variation=response::parseVariation(variationText);
  response::requireRunnableVariation(variation,emsContractText && *emsContractText);
  if(response::emsSelectionResponse(variation))ems_native::parseContract(emsContractText);
  const std::vector<std::string> paths = inputPaths(inputs);
  const Profile profile = configuredProfile(profileName, fileSystem(paths.front()), selection);
  const std::vector<CentralityClass> classes = centralityClasses(profile.system);
  const bool fused=dataJesUp || dataJesDown;
  std::vector<std::unique_ptr<DataJesView>> views;
  std::vector<DataJesView*> sharedViews;
  if(fused){
    if(!dataJesUp || !dataJesDown || std::string(dataJesUp)==dataJesDown || std::string(dataJesUp)==output || std::string(dataJesDown)==output)
      throw std::runtime_error("fused DATA JES requires three distinct output identities");
    if(profile.system!=System::pp || configuredDataJetScale!=1./.916 || !std::string(variationText).empty())
      throw std::runtime_error("fused DATA JES requires nominal pp DATA scale and empty response recipe");
    std::set<std::string> sources;
    for(const auto& path:paths){TFile source(path.c_str(),"READ");
      if(source.IsZombie() || fileSystem(source)!=System::pp || isSimulation(source) || !source.Get<TTree>("eventTree"))
        throw std::runtime_error("fused DATA JES requires original pp DATA eventTree inputs");
      if(!sources.insert(sourceSha(source)).second)throw std::runtime_error("duplicate original source in fused DATA JES");
    }
    for(double scale:{1./.907,1./.923}){views.push_back(std::make_unique<DataJesView>(scale,classes));sharedViews.push_back(views.back().get());}
  }
  Histograms h(classes);
  EventTable table, timeTable;
  const CentralityWeights weights = std::string(centralityWeightText).empty() ? CentralityWeights{} :
      centralityWeights(centralityWeightText, profile.system);
  std::unique_ptr<TFile> first(TFile::Open(paths.front().c_str(),"READ"));
  std::unique_ptr<Histograms> timeHists;
  if(!isSimulation(*first))timeHists=std::make_unique<Histograms>(classes);
  const bool retainTriggers=profile.configured && profile.system==System::auau && !isSimulation(*first);
  std::map<std::string,std::unique_ptr<Histograms>> triggerHists;
  std::map<std::string,EventTable> triggerTables;
  std::string nominalTrigger;
  std::unique_ptr<TFile> out(TFile::Open(output, "CREATE"));
  if (!out || out->IsZombie()) throw std::runtime_error(std::string("cannot create ") + output);
  std::unique_ptr<candidate_cache::Writer> cache;
  if(candidate_cache::enabled)cache=std::make_unique<candidate_cache::Writer>(*out,profile.jetRadius);
  if(retainTriggers) {
    nominalTrigger=SelectionJson::parse(profile.selectionJson).at("systems").at("auau").at("data_trigger").get<std::string>();
    for(const auto& key:{nominalTrigger}) triggerHists[key].reset(new Histograms(classes));
    Long64_t offset=0;
    for(const auto& key:{nominalTrigger}) for(size_t f=0;f<paths.size();++f) {
      fillFile(paths[f],static_cast<int>(f),profile,classes,*triggerHists.at(key),triggerTables[key],weights,key,offset,cache.get(),timeHists.get(),&timeTable);
      std::unique_ptr<TFile> source(TFile::Open(paths[f].c_str(),"READ"));
      offset+=static_cast<TTree*>(source->Get("eventTree"))->GetEntries();
    }
  } else {
    for (size_t f = 0; f < paths.size(); ++f) fillFile(paths[f], static_cast<int>(f), profile, classes, h, table, weights,"",0,cache.get(),timeHists.get(),&timeTable,&variation,emsContractText,digest,sharedViews);
  }
  out->cd();
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
  if(response::emsSelectionResponse(variation))TObjString(emsContractText).Write("ems_selection_contract");
  if(response::emsSelectionResponse(variation))TObjString(isSimulation(*first)?"true":"false").Write("ems_applied_to_simulation");
  if(timeHists){out->mkdir("timecut");out->GetDirectory("timecut")->cd();timeHists->write();timeTable.write();out->cd();}
  if(cache)cache->write();
  const auto nativeAncestry=photonjet_science::producerInputs(paths);
  photonjet_science::writeProducerEvidence(*out,nativeAncestry);
  if(fused){
    const char* targets[]={dataJesUp,dataJesDown};
    for(size_t i=0;i<views.size();++i){
      TFile varied(targets[i],"CREATE");if(varied.IsZombie())throw std::runtime_error("cannot create DATA JES output");
      views[i]->histograms.write();views[i]->records.write();
      TObjString(profile.name.c_str()).Write("profile");writeSelection(profile);
      varied.mkdir("timecut")->cd();views[i]->timing.write();views[i]->timingRecords.write();varied.cd();
      if(cache)for(auto* key: *out->GetListOfKeys()){
        const std::string name=key->GetName();if(name.rfind("table_candidates_",0)!=0)continue;
        auto* tree=out->Get<TTree>(name.c_str());varied.cd();auto* copy=tree->CloneTree(-1,"fast");copy->Write("",TObject::kOverwrite);
      }
      const SelectionJson proof={{"schema","PhotonJetDataJesOriginalTraversalV1"},{"view",i==0?"energy_up":"energy_down"},
        {"scale",views[i]->scale},{"science_readers",views[i]->scienceReaders},{"original_entries",views[i]->originalEntries},{"science_traversals",1},{"qa_extensions","nominal_reused"}};
      TObjString(proof.dump().c_str()).Write("data_jes_fusion");
      photonjet_science::writeProducerEvidence(varied,nativeAncestry);
      varied.Close();
    }
  }
  out->Close();
  return 0;
}
