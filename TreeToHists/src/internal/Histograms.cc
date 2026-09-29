// One adapter around the maintained physics fillers. Their macro entry points
// stay available for reference regressions and existing batch descriptors.
#include "Types.h"
#include "../FillPhotonJet.C"
#include "../FillPhotonQA.C"
#include "../FillPhotonJetResponse.C"
#include <THnSparse.h>
namespace photonjet {
static std::string numbers(const HistogramJob& values) {
  std::ostringstream out; out.precision(17);
  for (const auto& v: values) out << v.get<double>() << ' ';
  return out.str();
}
void fillHistogramJob(const HistogramJob& job) {
  const std::string in=job.at("input"), out=job.at("output"), sel=job.at("selection");
  const std::string weights=numbers(job.at("weights")), kind=job.at("kind");
  int code=0;
  if (kind=="spectra") code=FillPhotonJet(in.c_str(),out.c_str(),"nominal",sel.c_str(),weights.c_str());
  else if (kind=="qa") code=FillPhotonQA(in.c_str(),out.c_str(),"nominal",weights.c_str(),sel.c_str());
  else if (kind=="extensions") { fillExtensions(job); return; }
  else if (kind=="response") {
    const auto& h=job.at("histograms");
    code=FillPhotonJetResponse(in.c_str(),out.c_str(),"nominal",
      numbers(h.at("native_reco_response_edges_gev")).c_str(),
      numbers(h.at("native_truth_response_edges_gev")).c_str(),
      numbers(h.at("native_xj_edges")).c_str(),weights.c_str(),"",sel.c_str());
  } else throw std::runtime_error("unknown histogram operation: "+kind);
  if (code) throw std::runtime_error("ROOT filler failed: "+kind);
}

// Small additional QA pass. Missing optional stage-3 payloads are reported as
// unavailable, never represented by a histogram falsely claiming zero entries.
void fillExtensions(const HistogramJob& job) {
  std::unique_ptr<TFile> input(TFile::Open(job.at("input").get<std::string>().c_str()));
  if (!input || input->IsZombie()) throw std::runtime_error("extension input cannot open");
  const auto profile=validateSelectionInput(*input,configuredProfile("nominal",fileSystem(*input),job.at("selection").get<std::string>().c_str()));
  const bool auau=profile.system==System::auau, sim=isSimulation(*input);
  const DataEventSelection dataSelection(*input,profile);
  const auto classes=centralityClasses(profile.system);
  const auto weights=sim&&auau?centralityWeights(numbers(job.at("weights")),profile.system):CentralityWeights{};
  const auto schema=HistogramJob::parse(metaText(*input,"meta/schema"));
  std::vector<std::string> filled=schema.value("filled_reserved",std::vector<std::string>{});
  bool fractions=true;
  for (const auto* rr:{"02","03","04"}) for (const auto* det:{"emcal","ihcal","ohcal"})
    fractions=fractions&&std::find(filled.begin(),filled.end(),std::string("jet")+rr+"_"+det+"_frac")!=filled.end();
  if (job.at("histograms").value("require_stage3_fractions",false) && !fractions)
    throw std::runtime_error("stage-3 fraction payload required but not filled");
  qa::Book b(true);
  auto key=[](std::string x) { size_t p=0; while((p=x.find('/',p))!=std::string::npos) {x.replace(p,1,"__");p+=2;} return x; };
  const qa::Range cent=auau?qa::Range{100,0,100}:qa::Range{1,0,100};
  const qa::Range energy=auau?qa::Range{150,0,3000}:qa::Range{100,0,200};
  const qa::Range charge{150,0,3000}, mb{9,0,9}; // invalid, fail, pass for each of the two MB decisions
  auto h3=[&](std::string path,qa::Range x,qa::Range y,qa::Range z)->TH3F& {return b.h3(key(path),path.c_str(),x,y,z);};
  h3("qa/events/cutflow",{6,0,6},cent,{2,0,2});
  b.h2(key("qa/events/weights"),";event weight;finite/nonfinite",{200,-10,10},{2,0,2});
  b.h2(key("qa/events/mbd_charge"),";MBD charge;centrality",charge,cent);
  for(const auto* d:{"emcal","ihcal","ohcal","total"}) {
    h3(std::string("qa/events/")+d+"_vs_mbd_charge",energy,charge,mb);
    if(std::string(d)!="total") h3(std::string("qa/events/")+d+"_vs_centrality",energy,cent,mb);
  }
  for(const auto& c:classes) {
    if(fractions) for(const auto* rr:{"02","03","04"})
      h3("qa/jets/"+c.name+"/R"+rr+"/calo_fractions",{240,-1,2},{120,0,60},{3,0,3});
    if(sim) for(const auto* obj:{"photon","jet"})
      h3("qa/matching/"+c.name+"/"+obj,{100,0,1},{120,0,60},{8,0,8});
  }
  std::unique_ptr<THnSparseD> scoreRoles;
  if(sim) {
    int bins[]={100,30,auau?16:1,3}; double low[]={0,10,0,0}, high[]={1,40,auau?16.:1.,3};
    scoreRoles.reset(new THnSparseD(key("qa/photon_id/score_roles").c_str(),"score;ET;slice;unmatched,prompt,nonprompt",4,bins,low,high));
    scoreRoles->Sumw2();
  }
  TTreeReader r("eventTree",input.get());
  const VertexEventWeight vertexWeight(*input,r);
  TTreeReaderValue<int> run(r,"evt_run");
  TTreeReaderValue<ULong64_t> scaled(r,"trig_scaled");
  TTreeReaderValue<int> record(r,"evt_record"),ready(r,"evt_nodes_ready"),terminal(r,"evt_terminal_status"),valid(r,"vtx_valid"),cmb(r,"cent_mb_decision"),emb(r,"evt_embedded_mb_decision"),chargeState(r,"mbd_charge_sum_state");
  TTreeReaderValue<double> vz(r,"vtx_z"),c(r,"cent_value"),w(r,"w_event"),q(r,"mbd_charge_sum");
  std::vector<std::unique_ptr<TTreeReaderValue<double>>> calo;
  std::vector<std::unique_ptr<TTreeReaderValue<int>>> caloState;
  for(const auto* d:{"emcal","ihcal","ohcal","total"}) {
    calo.emplace_back(new TTreeReaderValue<double>(r,(std::string("calo_")+d).c_str()));
    caloState.emplace_back(new TTreeReaderValue<int>(r,(std::string("calo_")+d+"_state").c_str()));
  }
  std::vector<std::unique_ptr<TTreeReaderArray<double>>> pts,fracs;
  for(const auto* rr:{"02","03","04"}) {
    pts.emplace_back(new TTreeReaderArray<double>(r,(std::string("jet")+rr+"_pt").c_str()));
    if(fractions) for(const auto* d:{"emcal","ihcal","ohcal"})
      fracs.emplace_back(new TTreeReaderArray<double>(r,(std::string("jet")+rr+"_"+d+"_frac").c_str()));
  }
  TTreeReaderArray<double> et(r,"pho_et"),eta(r,"pho_eta"),phi(r,"pho_phi"),score(r,profile.scoreBranch.c_str());
  TTreeReaderArray<int> pre(r,"pho_presel_active"),state(r,profile.stateBranch.c_str());
  std::unique_ptr<TTreeReaderArray<int>> pr,pt,pc,jr,jc,tvalid,tclass;
  std::unique_ptr<TTreeReaderArray<double>> pm,jm;
  if(sim) {
    pr.reset(new TTreeReaderArray<int>(r,"plink_reco")); pt.reset(new TTreeReaderArray<int>(r,"plink_truth")); pc.reset(new TTreeReaderArray<int>(r,"plink_class"));
    jr.reset(new TTreeReaderArray<int>(r,"jlink_reco")); jc.reset(new TTreeReaderArray<int>(r,"jlink_class"));
    pm.reset(new TTreeReaderArray<double>(r,"plink_metric")); jm.reset(new TTreeReaderArray<double>(r,"jlink_metric"));
    tvalid.reset(new TTreeReaderArray<int>(r,"truthpho_g4_photon_valid"));tclass.reset(new TTreeReaderArray<int>(r,"truthpho_prompt_class"));
  }
  while(r.Next()) {
    reportProgress(r.GetCurrentEntry()+1);
    if(*record!=0) continue;
    if(!std::isfinite(*w)) throw std::runtime_error("nonfinite event weight");
    const double weight=*w*vertexWeight()*weights(*c), centValue=auau?*c:50.;
    auto& cut=h3("qa/events/cutflow",{6,0,6},cent,{2,0,2});
    const bool gates[]={true,*ready==1,*terminal==0,*valid==1,std::isfinite(*vz)&&std::abs(*vz)<profile.vertexAbsZMax,!auau||(*c>=0&&*c<80)};
    for(int i=0;i<6;++i) {if(!gates[i]) break;cut.Fill(i+.5,centValue,.5);cut.Fill(i+.5,centValue,1.5,weight);}
    b.h2(key("qa/events/weights"),"",{200,-10,10},{2,0,2}).Fill(*w,.5);
    if(!gates[1]||!gates[2]||!gates[3]) continue;
    const int mbState=(*cmb==1?2:*cmb==0?1:0)+3*(*emb==1?2:*emb==0?1:0);
    if(*chargeState==1&&std::isfinite(*q)) b.h2(key("qa/events/mbd_charge"),"",charge,cent).Fill(*q,centValue,weight);
    const char* names[]={"emcal","ihcal","ohcal","total"};
    for(int d=0;d<4;++d) if(**caloState[d]==1&&std::isfinite(**calo[d])) {
      if(*chargeState==1&&std::isfinite(*q)) h3(std::string("qa/events/")+names[d]+"_vs_mbd_charge",energy,charge,mb).Fill(**calo[d],*q,mbState+.5,weight);
      if(d<3) h3(std::string("qa/events/")+names[d]+"_vs_centrality",energy,cent,mb).Fill(**calo[d],centValue,mbState+.5,weight);
    }
    if(!gates[4]||!gates[5]||!dataSelection.accepts(*run,*scaled,*cmb)) continue;
    std::string cls="all";if(auau) for(const auto& cl:classes) if(*c>=cl.low&&*c<cl.high) cls=cl.name;
    if(fractions) for(int radius=0;radius<3;++radius) {
      const std::string rr=radius==0?"02":radius==1?"03":"04";
      for(int component=0;component<3;++component) {
        auto& a=*fracs[3*radius+component]; if(a.GetSize()!=pts[radius]->GetSize()) throw std::runtime_error("fraction alignment mismatch");
        for(size_t j=0;j<a.GetSize();++j) if(std::isfinite(a[j])&&std::isfinite((*pts[radius])[j]))
          h3("qa/jets/"+cls+"/R"+rr+"/calo_fractions",{240,-1,2},{120,0,60},{3,0,3}).Fill(a[j],(*pts[radius])[j],component+.5,weight);
      }
    }
    if(!sim) continue;
    std::vector<int> roles(et.GetSize(),0);
    for(size_t i=0;i<pr->GetSize();++i) {
      const int reco=(*pr)[i],truth=(*pt)[i]; if(reco<0) continue;
      if(size_t(reco)>=et.GetSize()) throw std::runtime_error("photon link index out of range");
      if(std::isfinite((*pm)[i])) h3("qa/matching/"+cls+"/photon",{100,0,1},{120,0,60},{8,0,8}).Fill((*pm)[i],et[reco],(*pc)[i]+.5,weight);
      if((*pc)[i]==0&&truth>=0) {
        if(size_t(truth)>=tclass->GetSize()) throw std::runtime_error("truth link out of range");
        if((*tvalid)[truth]==1) roles[reco]=((*tclass)[truth]==1||(*tclass)[truth]==2)?1:2;
      }
    }
    const int offset=pts[0]->GetSize()+pts[1]->GetSize();
    for(size_t i=0;i<jr->GetSize();++i) {
      const int j=(*jr)[i]-offset;
      if(j>=0&&size_t(j)<pts[2]->GetSize()&&std::isfinite((*jm)[i]))
        h3("qa/matching/"+cls+"/jet",{100,0,1},{120,0,60},{8,0,8}).Fill((*jm)[i],(*pts[2])[j],(*jc)[i]+.5,weight);
    }
    // The optional model's applicability is finite/NaN, independent of nominal state.
    for(size_t i=0;i<et.GetSize();++i) if(pre[i]==1&&(profile.optionalID||state[i]==1)&&std::isfinite(score[i])&&photonKinematics(profile,*vz,et[i],eta[i],phi[i])) {
      double x[]={score[i],et[i],auau?centralitySlice(*c)+.5:.5,roles[i]+.5}; scoreRoles->Fill(x,weight);
    }
  }
  if(r.GetEntryStatus()!=TTreeReader::kEntryBeyondEnd&&r.GetEntryStatus()!=TTreeReader::kEntryNotFound) throw std::runtime_error("extension read stopped before EOF");
  std::unique_ptr<TFile> out(TFile::Open(job.at("output").get<std::string>().c_str(),"CREATE"));
  if(!out||out->IsZombie()) throw std::runtime_error("cannot create extensions");
  b.write(); if(scoreRoles) scoreRoles->Write();
  TObjString(profile.name.c_str()).Write("profile");writeSelection(profile);
  HistogramJob availability={{"fractions",fractions?"available":"stage3_not_filled"},
    {"npb_timing","blocked_run_timing_calibration_and_tag_contract"},{"semantic_triggers","blocked_named_trigger_population_contract"}};
  TObjString(availability.dump().c_str()).Write("availability");out->Close();
}
}
