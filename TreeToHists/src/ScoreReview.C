#include "ScientificProducerEvidence.h"
/* ScoreReview: source-aware photon-ID score histograms from published eventTrees.
   Reuses the nominal native model reader and simulation ownership rules.
   The review uses the canonical newPPG12 preselection, one best cluster per matched
   prompt truth photon, and every selected inclusive-jet candidate. Training
   weights never enter this path. Source normalization is applied once when
   the complete, disjoint set of source products is collected. */
#include "PhotonJet.h"
#include "ScoreReviewScientific.h"
#include <TH1D.h>
#include <TSystem.h>
#include <set>
#include <iostream>

namespace photonjet_score_review {
using Json = nlohmann::json;
inline Json read(const std::string &path) {
  std::ifstream in(path); if (!in) throw std::runtime_error("cannot read score review " + path);
  Json value; in >> value; return value;
}
inline void require(bool test, const std::string &message) {
  if (!test) throw std::runtime_error(message);
}
inline void write(const std::string &path, const Json &value) {
  std::ofstream out(path+".partial"); out << value.dump(2) << '\n'; out.close();
  require(bool(out), "cannot write score review receipt");
  std::filesystem::rename(path+".partial",path);
}
inline int centralityClass(double cent) { return cent < 20 ? 0 : cent < 50 ? 1 : 2; }
inline int region(double score, double cent, const Json &cuts) {
  auto line=[&](int level) {const auto &p=cuts.at("WP"+std::to_string(level));
    return p.at("intercept").get<double>()+p.at("slope").get<double>()*cent;};
  return score > line(70) ? 3 : score > .1 && score < line(80) ? 1 : score > .1 ? 2 : 0;
}

int fill(const std::string &jobPath) {
  using namespace photonjet;
  const Json job=read(jobPath);
  require(job.at("schema")=="PhotonJetScoreReviewJobV1", "unsupported score review job");
  require(job.at("weighting")=="physics_event_times_sample_normalization", "score plots require physics weights");
  const std::string output=job.at("output");
  require(!std::filesystem::exists(output) && !std::filesystem::exists(output+".json"), "score review already exists");
  const auto &normalization=job.at("normalization");
  const bool photon=normalization.at("family")=="photonjet";
  require(photon || normalization.at("family")=="inclusivejet", "score review is simulation only");
  require(normalization.at("cross_section_pb").get<double>()>0, "invalid score review cross section");
  Ownership own;
  own.owner=photon ? Ownership::Owner::filterPhoton : Ownership::Owner::leadingTruthJet;
  own.low=normalization.at("window_gev")[0];
  if (!normalization.at("window_gev")[1].is_null()) own.high=normalization.at("window_gev")[1];
  Profile profile; profile.system=System::auau; profile.nativeID=true;
  profile.selectionJson=Json{{"systems",{{"auau",{{"optional_model",{{"native_evaluation",job.at("native_evaluation")}}}}}}}}.dump();
  const auto &cuts=job.at("working_points").at("cells").begin().value().at("cuts");
  std::vector<std::unique_ptr<TH1D>> histograms;
  for (int c=0;c<3;++c) for (const std::string label:{"signal","inclusive","truth_background"}) {
    auto h=std::make_unique<TH1D>(("score_"+std::to_string(c)+"_"+label).c_str(),";BDT score;event-weighted candidates",50,0,1);
    h->SetDirectory(nullptr); h->Sumw2();histograms.push_back(std::move(h));
  }
  Json counts=Json::array();
  for (int c=0;c<3;++c) counts.push_back({{"signal",0LL},{"inclusive",0LL},{"truth_background",0LL},
    {"regions",{{"below",0LL},{"non_tight",0LL},{"gap",0LL},{"tight",0LL}}}});
  const char *regionNames[]={"below","non_tight","gap","tight"};
  auto accumulate=[&](int c,int label,double value,double weight,double cent) {
    histograms.at(3*c+label)->Fill(std::min(value,std::nextafter(1.,0.)),weight);
    const char *names[]={"signal","inclusive","truth_background"};
    counts[c][names[label]]=counts[c][names[label]].get<long long>()+1;
    if (label<2) {const auto key=regionNames[region(value,cent,cuts)];
      counts[c]["regions"][key]=counts[c]["regions"][key].get<long long>()+1;}
  };
  progressPath=job.at("progress");progressIdentity=job.at("identity");progressTotal=job.at("expected_entries");
  progressOffset=0;progressStart=std::chrono::steady_clock::now();progressPhase="reading_score_histograms";
  progressRowTotal=0;entryBegin=0;reportProgress(0);
  Json sourceRecords=Json::array(),nativeAncestry=Json::array();std::set<std::string> sourceIds;
  long long denominator=0,ownedEvents=0;
  for (const auto &source:job.at("sources")) {
    const std::string path=source.at("path"),id=source.at("logical_id");
    require(sourceIds.insert(id).second,"duplicate score-review source identity");
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(),"READ"));
    require(file && !file->IsZombie(),"cannot open score-review source "+path);
    nativeAncestry.push_back(photonjet_science::producerInput(*file));
    require(fileSystem(*file)==System::auau && isSimulation(*file),"score review requires AuAu simulation");
    require(sourceSha(*file)==id,"score review logical source changed");
    const Json sourceMetadata=Json::parse(metaText(*file,"meta/source"));
    require(sourceMetadata.at("sample")==job.at("source_sample"),"score review sample differs");
    const Json weights=Json::parse(metaText(*file,"meta/weighting"));
    require(weights.value("composition","")=="unity" && !weights.value("absolute_normalization_applied",false) &&
      !weights.value("centrality_applied",false) && !weights.value("vertex_enabled",false),
      "historical before-reweighting review requires uncorrected producer event weights");
    auto *tree=file->Get<TTree>("eventTree");
    require(tree && tree->GetEntries()==source.at("entries").get<long long>(),"score-review entry count changed");
    TTreeReader reader(tree);OwnershipGate ownership(reader,own,true);
    TTreeReaderValue<int> record(reader,"evt_record"),ready(reader,"evt_nodes_ready"),valid(reader,"vtx_valid"),terminal(reader,"evt_terminal_status");
    TTreeReaderValue<double> cent(reader,"cent_value"),vz(reader,"vtx_z"),eventWeight(reader,"w_event");
    TTreeReaderArray<double> et(reader,"pho_et"),eta(reader,"pho_eta"),phi(reader,"pho_phi"),metric(reader,"plink_metric"),
      truthPt(reader,"truthpho_pt"),truthEta(reader,"truthpho_eta"),truthPhi(reader,"truthpho_phi");
    PhotonPreselection presel(reader);
    TTreeReaderArray<int> linkReco(reader,"plink_reco"),linkTruth(reader,"plink_truth"),
      linkClass(reader,"plink_class"),truthRole(reader,"truthpho_source_role"),barcode(reader,"truthpho_generator_barcode"),
      promptClass(reader,"truthpho_prompt_class"),dom(reader,"pho_dom_barcode"),domState(reader,"pho_dom_state"),domPid(reader,"pho_dom_pid");
    TTreeReaderArray<double> isoCone(reader,"isow_cone_sum"),isoThreshold(reader,"isow_threshold"),isoSideband(reader,"isow_sideband_threshold");
    TTreeReaderArray<int> isoMethod(reader,"isow_subtraction_method"),isoRole(reader,"isow_reconstructed_or_truth");
    PhotonScores score(*file,reader,profile,[&](size_t i) {return i<presel.GetSize() && presel[i]==1;});
    auto scientific=photonjet_science::qaContext(*file,profile,CentralityWeights(),DataEventSelection(*file,profile),0.,true,score.scientificFunctions());
    scientific["ownership"]={{"owner",photon?"filter_photon":"leading_truth_jet"},{"window",normalization.at("window_gev")}};
    for(int c=0;c<3;++c)for(int slot=0;slot<3;++slot) {
      const char* names[]={"signal","inclusive","truth_background"};auto& h=*histograms.at(3*c+slot);
      const auto incoming=photonjet_science::Trace::source(photonjet_science::scoreReviewDefinition(scientific,c,names[slot]));
      const auto old=photonjet_science::tag(&h,"physics_trace");
      photonjet_science::attachTrace(&h,old.is_null()?incoming:photonjet_science::Trace::combine("score_review_source_sum",{photonjet_science::Trace::fromJson(old),incoming}));
      photonjet_science::setTag(&h,"physics_covariance",photonjet_science::covariance("not_computed",{},0,"native event-weight Sumw2; cross-bin covariance not computed"));
    }
    progressRowTotal=tree->GetEntries();progressRowProcessed=0;
    long long localDenominator=0,localOwned=0,processed=0;
    while (reader.Next()) {
      ++processed;reportProgress(processed);
      const bool owned=ownership.owned();
      if (normalization.at("denominator")=="all_events" || owned) ++localDenominator;
      if (owned) ++localOwned;
      if (!owned || *record!=0 || *ready!=1 || *valid!=1 || *terminal!=0 ||
        !std::isfinite(*cent) || !(*cent>=0 && *cent<80) || !std::isfinite(*vz) || !(std::abs(*vz)<10)) continue;
      const double w=*eventWeight;
      require(std::isfinite(w) && w>=0,"invalid physics event weight");if (w==0) continue;
      const size_t n=et.GetSize(),nt=truthRole.GetSize();
      require(eta.GetSize()==n && phi.GetSize()==n && presel.GetSize()==n && score.GetSize()==n &&
        dom.GetSize()==n && domState.GetSize()==n && domPid.GetSize()==n && isoCone.GetSize()==2*n &&
        isoThreshold.GetSize()==2*n && isoSideband.GetSize()==2*n && isoMethod.GetSize()==2*n && isoRole.GetSize()==2*n,
        "score-review photon alignment differs");
      require(truthPt.GetSize()==nt && truthEta.GetSize()==nt && truthPhi.GetSize()==nt && barcode.GetSize()==nt &&
        promptClass.GetSize()==nt && linkReco.GetSize()==linkTruth.GetSize() && linkReco.GetSize()==metric.GetSize() &&
        linkReco.GetSize()==linkClass.GetSize(),"score-review truth/link alignment differs");
      std::vector<int> matched(n,-1);std::vector<double> bestMetric(n,std::numeric_limits<double>::infinity());
      for (size_t l=0;l<linkReco.GetSize();++l) {
        if (linkClass[l]!=kLinkMatch) continue;
        const int r=linkReco[l],t=linkTruth[l];require(r>=0 && size_t(r)<n && t>=0 && size_t(t)<nt,"invalid prompt-match indices");
        if (truthRole[t]!=kSignalSourceRole || !std::isfinite(truthPt[t]) || !(truthPt[t]>0) ||
          !std::isfinite(truthEta[t]) || !(std::abs(truthEta[t])<kTruthPhotonAbsEtaMax) || !std::isfinite(truthPhi[t])) continue;
        require(std::isfinite(metric[l]),"nonfinite final prompt-match metric");
        require(matched[r]<0 || matched[r]==t,"one candidate has two final prompt truth matches");
        if (metric[l]<bestMetric[r]) {matched[r]=t;bestMetric[r]=metric[l];}
      }
      std::map<int,int> winners;
      std::map<int,std::set<int>> truthClasses;
      for (size_t t=0;t<nt;++t) truthClasses[barcode[t]].insert(promptClass[t]);
      const int c=centralityClass(*cent);
      for (size_t i=0;i<n;++i) {
        if (presel[i]!=1 || !std::isfinite(et[i]) || !(et[i]>=15 && et[i]<35) || !std::isfinite(eta[i]) ||
          !(std::abs(eta[i])<.7) || !std::isfinite(phi[i])) continue;
        // The historical candidate cache required a finite reconstructed SUB1
        // R=0.4 witness, without selecting isolated photons.
        const size_t iso=2*i+1;
        if (isoMethod[iso]!=2 || isoRole[iso]!=0 || !std::isfinite(isoCone[iso]) ||
          !std::isfinite(isoThreshold[iso]) || !std::isfinite(isoSideband[iso])) continue;
        const double s=score[i];if (!std::isfinite(s)) continue;
        require(s>=0 && s<=1,"invalid saved-model score");
        if (photon) {
          if (matched[i]<0) continue;
          auto found=winners.find(matched[i]);
          if (found==winners.end() || std::make_tuple(bestMetric[i],-et[i],i)<
            std::make_tuple(bestMetric[found->second],-et[found->second],size_t(found->second))) winners[matched[i]]=i;
        } else {
          accumulate(c,1,s,w,*cent);
          const auto &classes=truthClasses[dom[i]];
          if (domState[i]==2 && classes.size()<=1 && !(domPid[i]==22 && classes.empty()) &&
            (classes.empty() || *classes.begin()>=1) && !(classes.size()==1 && (*classes.begin()==1 || *classes.begin()==2)))
            accumulate(c,2,s,w,*cent);
        }
      }
      for (const auto &[truth,i]:winners) accumulate(c,0,score[i],w,*cent);
    }
    require(processed==tree->GetEntries(),"score review did not read its complete source");
    denominator+=localDenominator;ownedEvents+=localOwned;
    sourceRecords.push_back(Json{{"logical_id",id},{"sha256",source.at("sha256")},{"entries",processed},
      {"denominator_events",localDenominator},{"owned_events",localOwned}});
    progressOffset+=processed;
  }
  require(progressOffset==progressTotal,"score-review work differs from pinned total");
  TFile product((output+".partial").c_str(),"CREATE");require(!product.IsZombie(),"cannot create score-review product");
  for (auto &h:histograms) h->Write();
  Json receipt={{"schema","PhotonJetScoreReviewProductV1"},{"status","PASS_EVENT_WEIGHTED_PENDING_SAMPLE_NORMALIZATION"},
    {"identity",job.at("identity")},{"weighting",job.at("weighting")},{"source_sample",job.at("source_sample")},
    {"model",job.at("model")},{"normalization",normalization},{"denominator_events",denominator},{"owned_events",ownedEvents},
    {"classes",counts},{"sources",sourceRecords},{"working_points",job.at("working_points")},
    {"selection",{{"preselection",canonicalPreselection()},{"signal","one best cluster per event and source-role-1 matched truth photon"},
      {"inclusive","every selected inclusive-jet candidate, without truth-class filtering"},{"et_gev",{15,35}},
      {"eta_abs_max",.7},{"vertex_abs_max_cm",10},{"centrality",{0,80}},
      {"winner_order",{"smallest plink_metric","largest pho_et","lowest candidate index"}},
      {"isolation_cut",false},{"finite_reconstructed_sub1_r04_witness",true},
      {"centrality_reweighting",false},{"vertex_reweighting",false},{"training_weights_used",false}}}};
  TObjString(receipt.dump().c_str()).Write("review");
  std::vector<photonjet_science::Trace> reviewInputs;for(const auto& h:histograms)reviewInputs.push_back(photonjet_science::traceOf(h.get()));
  TObjString(Json{{"trace",photonjet_science::Trace::combine("score_review_region_counts",reviewInputs,{{"cuts",cuts}}).json()}}.dump().c_str()).Write("physics_review");
  photonjet_science::writeProducerEvidence(product,nativeAncestry);
  product.Close();
  std::filesystem::rename(output+".partial",output);write(output+".json",receipt);
  progressPhase="complete";progressRowTotal=0;reportProgress(0);
  std::cout << "[ScoreReview] complete " << progressTotal << " entries, " << sourceRecords.size() << " sources\n";
  return 0;
}

/* Collection verifies exact source coverage before applying each sample's
   ownership denominator. Neither ROOT hadd nor a histogram-only sum proves
   that the inputs are disjoint or complete. */
int collect(const std::string &jobPath) {
  const Json job=read(jobPath);require(job.at("schema")=="PhotonJetScoreReviewCollectionJobV1","unsupported score collection");
  require(job.at("weighting")=="physics_event_times_sample_normalization","score plots require physics weights");
  const std::string output=job.at("output");
  require(!std::filesystem::exists(output) && !std::filesystem::exists(output+".json"),"score collection already exists");
  photonjet::progressPath=job.value("progress","");photonjet::progressIdentity=job.value("identity",Json::object());
  photonjet::progressTotal=job.at("products").size();photonjet::progressRowTotal=photonjet::progressTotal;
  photonjet::progressOffset=0;photonjet::entryBegin=0;photonjet::progressStart=std::chrono::steady_clock::now();
  photonjet::progressPhase="collecting_score_histograms";photonjet::reportProgress(0);
  std::map<std::string,Json> expected;
  for (const auto &s:job.at("sources")) require(expected.emplace(s.at("logical_id").get<std::string>(),s).second,"duplicate expected score source");
  std::set<std::string> observed;
  std::map<std::string,std::vector<std::unique_ptr<TH1D>>> samples;
  std::map<std::string,long long> denominators;
  std::map<std::string,Json> normalizations;
  Json counts=Json::array(),sourceRecords=Json::array(),selection;long long completed=0;
  for (int c=0;c<3;++c) counts.push_back({{"low",c==0 ? 0 : c==1 ? 20 : 50},{"high",c==0 ? 20 : c==1 ? 50 : 80},
    {"signal",0LL},{"background",0LL},{"truth_background",0LL},{"rows",0LL},
    {"regions",{{"below",0LL},{"non_tight",0LL},{"gap",0LL},{"tight",0LL}}}});
  Json nativeAncestry=Json::array();
  for (const auto &input:job.at("products")) {
    std::unique_ptr<TFile> f(TFile::Open(input.get<std::string>().c_str(),"READ"));require(f && !f->IsZombie(),"missing score product");
    nativeAncestry.push_back(photonjet_science::producerInput(*f));
    const auto *meta=dynamic_cast<TObjString*>(f->Get("review"));require(meta,"missing score-product receipt");
    const Json p=Json::parse(meta->GetString().Data());const std::string sample=p.at("source_sample");
    require(p.at("schema")=="PhotonJetScoreReviewProductV1" && p.at("model")==job.at("model") &&
      p.at("working_points")==job.at("working_points") && p.at("weighting")==job.at("weighting"),"score-product model/cuts/weights differ");
    if (selection.is_null()) selection=p.at("selection");
    require(selection==p.at("selection"),"score-product candidate selection differs");
    if (normalizations.count(sample)) require(normalizations[sample]==p.at("normalization"),"source normalization differs");
    normalizations[sample]=p.at("normalization");denominators[sample]+=p.at("denominator_events").get<long long>();
    for (const auto &s:p.at("sources")) {
      const std::string id=s.at("logical_id");require(observed.insert(id).second,"duplicate source at score merge");
      require(expected.count(id) && expected.at(id).at("sha256")==s.at("sha256") && expected.at(id).at("entries")==s.at("entries") &&
        expected.at(id).at("sample")==sample,"score source differs from exact input snapshot");sourceRecords.push_back(s);
    }
    auto &sum=samples[sample];int slot=0;
    for (int c=0;c<3;++c) for (const std::string label:{"signal","inclusive","truth_background"}) {
      auto *h=dynamic_cast<TH1D*>(f->Get(("score_"+std::to_string(c)+"_"+label).c_str()));require(h && h->GetSumw2N()>0,"score product lacks histogram uncertainties");
      if (sum.size()==size_t(slot)) {sum.emplace_back(static_cast<TH1D*>(h->Clone()));sum.back()->SetDirectory(nullptr);}
      else {const auto trace=photonjet_science::Trace::combine("score_review_product_sum",{photonjet_science::traceOf(sum[slot].get()),photonjet_science::traceOf(h)});
        require(sum[slot]->Add(h),"incompatible score binning");photonjet_science::attachTrace(sum[slot].get(),trace);}++slot;
      const std::string key=label=="inclusive" ? "background" : label;
      counts[c][key]=counts[c][key].get<long long>()+p.at("classes")[c].at(label).get<long long>();
    }
    for (int c=0;c<3;++c) for (const std::string key:{"below","non_tight","gap","tight"})
      counts[c]["regions"][key]=counts[c]["regions"][key].get<long long>()+p.at("classes")[c].at("regions").at(key).get<long long>();
    photonjet::reportProgress(++completed);
  }
  require(observed.size()==expected.size(),"incomplete score source coverage");
  std::vector<std::unique_ptr<TH1D>> total;Json factors=Json::object();
  for (auto &[sample,hs]:samples) {
    require(denominators[sample]>0,"empty simulation normalization denominator");
    const double factor=normalizations[sample].at("cross_section_pb").get<double>()/denominators[sample];
    factors[sample]={{"cross_section_pb",normalizations[sample].at("cross_section_pb")},{"denominator_events",denominators[sample]},
      {"pb_per_event",factor},{"definition",normalizations[sample]}};
    for (size_t slot=0;slot<hs.size();++slot) {
      hs[slot]->Scale(factor);
      photonjet_science::attachTrace(hs[slot].get(),photonjet_science::scoreReviewNormalization(
          photonjet_science::traceOf(hs[slot].get()),sample,normalizations[sample].at("cross_section_pb"),denominators[sample]));
      if (total.size()==slot) {total.emplace_back(static_cast<TH1D*>(hs[slot]->Clone()));total.back()->SetDirectory(nullptr);}
      else {const auto trace=photonjet_science::Trace::combine("score_review_sample_sum",{photonjet_science::traceOf(total[slot].get()),photonjet_science::traceOf(hs[slot].get())});
        require(total[slot]->Add(hs[slot].get()),"incompatible score sample bins");photonjet_science::attachTrace(total[slot].get(),trace);}
    }
  }
  long long rows=0;for (auto &c:counts) {c["rows"]=c["signal"].get<long long>()+c["background"].get<long long>();rows+=c["rows"].get<long long>();}
  Json info={{"schema","PhotonIDWorkingPointReviewV1"},{"status","PASS"},{"rows",rows},{"classes",counts},
    {"weighting",job.at("weighting")},{"working_points",job.at("working_points")},{"sample_normalization",factors},
    {"population","historically preselected exclusive prompt winners and all inclusive-jet candidates from exact published simulation"},
    {"provenance",job.at("provenance")},{"sources",sourceRecords},{"selection",selection},{"training_weights_used",false}};
  TFile out((output+".partial").c_str(),"CREATE");require(!out.IsZombie(),"cannot create score collection");
  for (int c=0;c<3;++c) {
    total[3*c]->Write(("score_"+std::to_string(c)+"_signal").c_str());
    total[3*c+1]->Write(("score_"+std::to_string(c)+"_background").c_str());
    total[3*c+2]->Write(("score_"+std::to_string(c)+"_truth_background").c_str());
  }
  TObjString(info.dump().c_str()).Write("review");
  std::vector<photonjet_science::Trace> reviewInputs;for(const auto& h:total)reviewInputs.push_back(photonjet_science::traceOf(h.get()));
  TObjString(Json{{"trace",photonjet_science::Trace::combine("collected_score_review_region_counts",reviewInputs,
      {{"cuts",job.at("working_points").at("cells").begin().value().at("cuts")}}).json()}}.dump().c_str()).Write("physics_review");
  photonjet_science::writeProducerEvidence(out,nativeAncestry);
  out.Close();std::filesystem::rename(output+".partial",output);write(output+".json",info);
  photonjet::progressPhase="complete";photonjet::reportProgress(completed);return 0;
}
} // namespace photonjet_score_review

void ScoreReview(const char *job) {
  try {const auto spec=photonjet_score_review::read(job);
    if (spec.at("schema")=="PhotonJetScoreReviewCollectionJobV1") photonjet_score_review::collect(job);
    else photonjet_score_review::fill(job);}
  catch (const std::exception &e) {std::cerr << "ScoreReview: " << e.what() << '\n';gSystem->Exit(1);}
}
