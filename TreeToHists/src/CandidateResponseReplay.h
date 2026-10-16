#pragma once
/* Native response replay. The caller binds source identities, current
   selection and complete sample normalization to the accepted package. */
#include "CandidateMeasurementReplay.h"
#include "FillPhotonJetResponse.C"
#include "CentralityReplay.h"

namespace photonjet { namespace candidate_replay {
inline Summary responseMeasurement(TFile& input,const std::string& prefix,const Profile& profile,
    const std::vector<CentralityClass>& classes,const response::ResponseGrids& grids,
    std::vector<std::vector<response::Accumulators>>& output,const std::vector<std::string>& sourceIdentities,
    double normalization,const response::Variation& variation,Long64_t maxCandidates=2000000,
    const energy_selection::Context* energySelection=nullptr,const Slice* slice=nullptr,
    response::ResolutionMigration* migration=nullptr, const nlohmann::json* admittedCentralityWeights=nullptr) {
  response::validatePhotonSelection(variation);
  response::requireRunnableVariation(variation);
  const bool only=response::responseOnly(variation);
  require(!only||!energySelection,"response-only resolution must retain original classification");
  const bool photonChanged=variation.photonScale!=1. || variation.emScale!=1. || variation.photonSmear!=0. || variation.photonResolution!="off";
  require(!photonChanged||energySelection||only,"photon response replay variation requires whole-selection evaluator");
  require(profile.simulation,"response replay requires simulation");
  require(!variation.centralityOff || (profile.system == System::auau && admittedCentralityWeights),
          "centrality off requires admitted AuAu weights");
  const int radius=candidate_cache::radiusIndex(profile.jetRadius);
  require(prefix=="records/candidates/" || prefix=="table_candidates_","unsupported response table layout");
  require(!sourceIdentities.empty() && std::isfinite(normalization) && normalization>0 && maxCandidates>0,
          "invalid response replay envelope");
  require(profile.system!=System::auau || (profile.isolationRadiusIndex==0 && profile.isolationSubtraction==2),
          "AuAu response requires R03 SUB1 isolation");
  require(profile.system!=System::pp || (profile.isolationRadiusIndex==1 && profile.isolationSubtraction==1),
          "pp response requires R04 topocluster isolation");
  candidate_replay::sourceIdentities(sourceIdentities,slice);
  const bool packaged=prefix=="records/candidates/";
  auto* events=input.Get<TTree>((prefix+"events").c_str());require(events,"missing response event table");
  for(const auto& context:photonjet_science::replayScientificContexts(*events,profile,slice?slice->source:-1,
      nlohmann::json::parse(response::variationJson(variation)),{{"reco",grids.reco},{"truth",grids.truth},{"xj",grids.xj}},
      normalization,false,admittedCentralityWeights,variation.centralityOff))
    for(auto& group:output)for(auto& accumulator:group)accumulator.scientificContext(context);
  events->SetBranchStatus("*",false);
  std::map<std::string,TLeaf*> ei,ed;
  for(const char* name:{"entry","file","photon_input_count","photon_retained_count","photon_below_floor_count","photon_invalid_et_count"})
    ei[name]=leaf(*events,name,"Long64_t");
  for(const char* kind:{"jet","jetunsub","truthjet"})for(const char* radius:{"02","03","04"}) {
    const std::string name=std::string(kind)+radius+"_count";ei[name]=leaf(*events,name,"Long64_t");
  }
  if(packaged)ei["source"]=leaf(*events,"source","Long64_t");
  ei["response_match_witness"]=leaf(*events,"response_match_witness","Long64_t");
  for(const char* name:{"centrality","vtx_z","weight","weight_base","weight_vertex","weight_centrality","reco_et_min_gev","jet_radius"})
    ed[name]=leaf(*events,name,"Double_t");
  std::vector<std::string> pf={"pho_et","pho_eta","pho_phi","selected_score","applied_preselection"};
  if(!profile.optionalID)pf.push_back(profile.stateBranch);
  if(profile.augmentedClasses){pf.push_back("pho_id_class");pf.push_back("pho_iso_class");}
  pf.push_back("response_truth_index");pf.push_back("response_match_metric");
  Rows photons(input,prefix,"photons",pf,packaged,energySelection?energy_selection::capturedFields(profile):std::vector<std::string>{},slice);
  Rows isolation(input,prefix,"isolation",{"isow_cone_sum","isow_threshold","isow_sideband_threshold",
      "isow_subtraction_method","isow_reconstructed_or_truth"},packaged,{},slice);
  std::vector<std::string> jf={"jet_pt","jet_eta","jet_phi","jet_rawpt"};
  if(variation.emScale!=1.)jf.push_back("jet_emcal_frac");
  jf.push_back("response_truth_index");
  Rows jets(input,prefix,"jets",jf,packaged,{},slice);
  Rows truths(input,prefix,"truth_photons",{"truthpho_pt","truthpho_eta","truthpho_phi","truthpho_source_role",
      "truthpho_prompt_class","truthpho_g4_photon_valid","truthpho_truth_isolation_r03"},packaged,{},slice);
  Rows truthJets(input,prefix,"truth_jets",{"truthjet_pt","truthjet_eta","truthjet_phi"},packaged,{},slice);
  for(const char* name:{"unsub_jets","photon_links","jet_links"}){
    auto* table=input.Get<TTree>((prefix+name).c_str());require(table && table->GetEntries()==0,"mixed raw links and verified match representation");
  }
  Summary summary;Long64_t previousSource=-1,previousEntry=-1,previousFile=-1;
  const double nan=std::numeric_limits<double>::quiet_NaN();
  const auto eventRange=rangeOf(*events,"events",slice);
  for(Long64_t event=eventRange.first;event<eventRange.second;++event) {
    require(events->GetEntry(event)>0,"unreadable response event");
    const auto i=[&](const std::string& key){return ei.at(key)->GetValueLong64();};
    const auto d=[&](const std::string& key){return ed.at(key)->GetValue();};
    require(i("response_match_witness")==1,"SIM record lacks its verified match witness");
    require(d("jet_radius")==profile.jetRadius,"record jet radius differs from response profile");
    const Long64_t source=packaged?i("source"):i("file"),entry=i("entry"),file=i("file");
    require(source>=0 && source<static_cast<Long64_t>(sourceIdentities.size()) && source>=previousSource && entry>=0 && file>=0,
            "invalid original response source/entry");
    require(!slice || source==slice->source,"response event outside requested source slice");
    if(source==previousSource && file==previousFile)require(entry>previousEntry,"duplicate response event");
    previousSource=source;previousEntry=entry;previousFile=file;
    const auto total=i("photon_input_count"),retained=i("photon_retained_count"),dropped=i("photon_below_floor_count"),invalid=i("photon_invalid_et_count");
    require(total>=0 && total<=maxCandidates && retained>=0 && dropped>=0 && invalid>=0 && total==retained+dropped+invalid,
            "response photon retention accounting differs");
    const double floor=d("reco_et_min_gev"),scale=variation.photonScale*variation.emScale;
    require(std::isfinite(floor) && floor>=0 && scale>0 && std::isfinite(scale),"invalid photon replay floor/scale");
    require(!dropped || (variation.photonSmear==0 && variation.photonResolution=="off" && floor*scale<=profile.photonEtMin),
            "discarded photons could migrate into the varied selection; retain all finite simulation photons");
    const double weight=d("weight"),product=d("weight_base")*d("weight_vertex")*d("weight_centrality");
    require(std::isfinite(weight) && std::isfinite(product) && weight==product && std::isfinite(weight*normalization),
            "response event-weight factors differ");
    if (admittedCentralityWeights)
      photonjet_centrality_replay::recordedFactor(d("centrality"), d("weight_centrality"), *admittedCentralityWeights);
    const double responseWeight = variation.centralityOff ?
      centrality_source::recordedWeight(weight, d("weight_base"), d("weight_vertex"), d("weight_centrality"), normalization, true) : weight * normalization;
    require(responseEvent(profile,d("vtx_z")),"response event outside original vertex acceptance");
    auto ps=photons.take(event,source,maxCandidates),qs=isolation.take(event,source,2*maxCandidates),
         js=jets.take(event,source,maxCandidates),ts=truths.take(event,source,maxCandidates),tjs=truthJets.take(event,source,maxCandidates);
    require(ps.size()==static_cast<size_t>(retained) && qs.size()==2*ps.size(),"response photon/isolation cardinality differs");
    const auto counts=[&](const char* kind){std::array<size_t,3> result{};int r=0;
      for(const char* radius:{"02","03","04"}){const auto n=i(std::string(kind)+radius+"_count");
        require(n>=0 && n<=maxCandidates,"invalid jet collection count");result[r++]=n;}return result;};
    const auto ranges=response::JetRanges::of(counts("jet"),counts("jetunsub"),counts("truthjet"),radius);
    require(js.size()==ranges.nJet && tjs.size()==ranges.nTruthJet,
            "response jet collection cardinality differs");
    Values<double> stored(total,nan),eta(total,nan),phi(total,nan),score(total,nan),
      cone(2*total,nan),threshold(2*total,nan),sideband(2*total,nan);
    Values<int> presel(total,0),state(total,0),id(total,0),iso(total,0),method(2*total,0),role(2*total,0);
    for(size_t k=0;k<ps.size();++k) {
      const auto& p=ps[k];const auto original=integer(p.at("__input_index"),"photon input index");
      require(original>=0 && original<total,"response photon index outside original collection");
      require(qs[2*k].at("__input_index")==2*original && qs[2*k+1].at("__input_index")==2*original+1,
              "response isolation no longer belongs to photon");
      require(std::isfinite(p.at("pho_et")) && p.at("pho_et")>=floor,"response photon below recorded floor");
      stored[original]=p.at("pho_et");eta[original]=p.at("pho_eta");phi[original]=p.at("pho_phi");score[original]=p.at("selected_score");
      presel[original]=integer(p.at("applied_preselection"),"preselection");
      require(presel[original]==0 || presel[original]==1,"invalid response preselection");
      if(!profile.optionalID)state[original]=integer(p.at(profile.stateBranch),"score state");
      if(profile.augmentedClasses){id[original]=integer(p.at("pho_id_class"),"ID class");iso[original]=integer(p.at("pho_iso_class"),"isolation class");}
      for(size_t radius=0;radius<2;++radius){const auto& q=qs[2*k+radius];const size_t at=2*original+radius;
        cone[at]=q.at("isow_cone_sum");threshold[at]=q.at("isow_threshold");sideband[at]=q.at("isow_sideband_threshold");
        method[at]=integer(q.at("isow_subtraction_method"),"isolation method");role[at]=integer(q.at("isow_reconstructed_or_truth"),"isolation role");}
    }
    Values<double> jetPt,jetEta,jetPhi,rawPt,em,truthPt,truthEta,truthPhi,truthIso,truthJetPt,truthJetEta,truthJetPhi;
    Values<int> truthRole,truthPrompt,truthValid;
    const auto identity=[](const auto& rows){for(size_t k=0;k<rows.size();++k)
      require(rows[k].at("__input_index")==static_cast<double>(k),"unfiltered response collection order changed");};
    identity(js);identity(ts);identity(tjs);
    for(const auto& j:js){jetPt.push_back(j.at("jet_pt"));jetEta.push_back(j.at("jet_eta"));jetPhi.push_back(j.at("jet_phi"));rawPt.push_back(j.at("jet_rawpt"));
      if(variation.emScale!=1.)em.push_back(j.at("jet_emcal_frac"));}
    for(const auto& t:ts){truthPt.push_back(t.at("truthpho_pt"));truthEta.push_back(t.at("truthpho_eta"));truthPhi.push_back(t.at("truthpho_phi"));
      truthRole.push_back(integer(t.at("truthpho_source_role"),"truth role"));truthPrompt.push_back(integer(t.at("truthpho_prompt_class"),"truth prompt"));
      truthValid.push_back(integer(t.at("truthpho_g4_photon_valid"),"truth valid"));truthIso.push_back(t.at("truthpho_truth_isolation_r03"));}
    for(const auto& t:tjs){truthJetPt.push_back(t.at("truthjet_pt"));truthJetEta.push_back(t.at("truthjet_eta"));truthJetPhi.push_back(t.at("truthjet_phi"));}
    response::EventMatches match;
    {
      match={std::vector<long>(total,-1),std::vector<double>(total,nan),std::vector<long>(ranges.nJet,-1)};
      for(const auto& p:ps){
        const auto original=integer(p.at("__input_index"),"original photon match index"),truth=integer(p.at("response_truth_index"),"photon truth match");
        require(truth>=-1 && truth<static_cast<Long64_t>(ts.size()),"verified photon match outside truth collection");
        match.photonTruth.at(original)=truth;match.photonMetric.at(original)=p.at("response_match_metric");
      }
      std::set<long> usedTruth;
      for(size_t k=0;k<js.size();++k){
        const auto truth=integer(js[k].at("response_truth_index"),"jet truth match");
        require(truth>=-1 && truth<static_cast<Long64_t>(tjs.size()),"verified jet match outside truth collection");
        require(truth<0 || usedTruth.insert(truth).second,"duplicate verified jet truth match");match.jetTruth[k]=truth;
      }
    }
    const auto key=response::textKey(sourceIdentities[source]);
    response::Varied et(stored,variation.photonResolution=="off"?scale:1.,variation.photonSmear,variation.seed,key,entry,0);
    response::Varied variedJets(jetPt,variation.jetResolution=="off"?variation.jetScale:1.,variation.jetSmear,variation.seed,key,entry,1);
    if(variation.photonResolution!="off"){et.smearPhotons(truthPt,match.photonTruth,variation,key,entry);et.scaleAll(scale);}
    if(variation.jetResolution!="off"){variedJets.smearJets(truthJetPt,match.jetTruth,variation,profile.jetRadius,key,entry);variedJets.scaleAll(variation.jetScale);}
    if(variation.emScale!=1.)variedJets.scaleEmFraction(em,variation.emScale);
    if(energySelection)for(size_t k=0;k<ps.size();++k) {
      const auto original=integer(ps[k].at("__input_index"),"photon input index");
      auto input=energy_selection::recordInput(ps[k],{qs[2*k],qs[2*k+1]},d("vtx_z"),d("centrality"));
      const auto result=energySelection->apply(input,et[original]);
      score[original]=result.score;state[original]=result.scoreState;presel[original]=result.preselection;
      for(size_t r=0;r<2;++r){cone[2*original+r]=input.prescribedCone[r];threshold[2*original+r]=result.threshold[r];sideband[2*original+r]=result.sideband[r];}
    }
    for(size_t k=0;k<ts.size();++k)
      require(analysisTruthSignal(truthValid[k],truthPrompt[k],truthEta[k],truthIso[k])==(truthRole[k]==kSignalSourceRole),
              "cached truth signal definition differs");
    response::accumulateResponse(profile,classes,grids,output,d("vtx_z"),d("centrality"),responseWeight,
        et,eta,phi,score,presel,state,cone,threshold,sideband,method,role,id,iso,variedJets,jetEta,jetPhi,
        truthPt,truthEta,truthPhi,truthRole,truthJetPt,truthJetEta,truthJetPhi,match,
        only?static_cast<const std::vector<double>*>(&stored):nullptr,migration);
    ++summary.cacheEvents;reportProgress(summary.cacheEvents);
  }
  photons.finish();isolation.finish();jets.finish();truths.finish();truthJets.finish();
  return summary;
}
} }
