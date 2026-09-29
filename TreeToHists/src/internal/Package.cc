// Consolidate source products without losing source bootstrap contributions or
// event references. At most one input product and one accumulated sample are
// resident; contribution and event tables stream directly to ROOT baskets.
#include "Types.h"
#include "../EventTreeInput.h"
#include <TKey.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TLeaf.h>
#include <TTree.h>
#include <map>
#include <memory>

namespace photonjet {
namespace packaging {
TDirectory* directory(TFile& file,const std::string& path) {
  if(path.empty()) return &file;
  if(!file.GetDirectory(path.c_str())) file.mkdir(path.c_str(),"",true);
  auto* d=file.GetDirectory(path.c_str());if(!d) throw std::runtime_error("cannot create package directory");return d;
}
void sameAxes(const TH1& a,const TH1& b) {
  if(a.GetDimension()!=b.GetDimension()||a.GetNcells()!=b.GetNcells()) throw std::runtime_error("histogram shape mismatch");
  const TAxis* aa[]={a.GetXaxis(),a.GetYaxis(),a.GetZaxis()};
  const TAxis* bb[]={b.GetXaxis(),b.GetYaxis(),b.GetZaxis()};
  for(int d=0;d<a.GetDimension();++d) {
    if(aa[d]->GetNbins()!=bb[d]->GetNbins()) throw std::runtime_error("histogram axis mismatch");
    for(int i=1;i<=aa[d]->GetNbins()+1;++i) if(aa[d]->GetBinLowEdge(i)!=bb[d]->GetBinLowEdge(i)) throw std::runtime_error("histogram edges mismatch");
  }
}
struct Table {
  std::unique_ptr<TTree> output;
  std::map<std::string,Long64_t> ints;
  std::map<std::string,Double_t> doubles;
  void append(TTree& input,const std::string& name,Long64_t source,Long64_t offset,TDirectory* dest) {
    if(!output) {
      dest->cd();output.reset(new TTree(name.c_str(),input.GetTitle()));output->SetAutoFlush(-8*1024*1024);
      for(auto* obj:*input.GetListOfLeaves()) {
        auto* leaf=static_cast<TLeaf*>(obj);const std::string n=leaf->GetName(),t=leaf->GetTypeName();
        if(t=="Double_t"||t=="Float_t") {doubles[n]=0;output->Branch(n.c_str(),&doubles[n],(n+"/D").c_str());}
        else {ints[n]=0;output->Branch(n.c_str(),&ints[n],(n+"/L").c_str());}
      }
      ints["source"]=0;output->Branch("source",&ints["source"],"source/L");
    }
    for(Long64_t row=0;row<input.GetEntries();++row) {
      if(input.GetEntry(row)<=0) throw std::runtime_error("cannot read event contribution");
      for(auto& item:doubles) {auto* l=input.GetLeaf(item.first.c_str());if(!l) throw std::runtime_error("table schema changed");item.second=l->GetValue();}
      for(auto& item:ints) {
        if(item.first=="source") {item.second=source;continue;}
        auto* l=input.GetLeaf(item.first.c_str());if(!l) throw std::runtime_error("table schema changed");
        item.second=l->GetValueLong64();
        if(item.first=="event") item.second+=offset;
        if(item.first=="file") item.second=source;
      }
      output->Fill();
    }
  }
  void write() {if(output) {output->GetDirectory()->cd();output->Write("",TObject::kOverwrite);output->SetDirectory(nullptr);}}
};
}
void packageHistogramJob(const HistogramJob& job) {
  using namespace packaging;
  std::unique_ptr<TFile> out(TFile::Open(job.at("output").get<std::string>().c_str(),"CREATE"));
  if(!out||out->IsZombie()) throw std::runtime_error("cannot create sample package");
  auto* records=directory(*out,"records");records->cd();
  progressPath=job.value("progress","");progressTotal=job.at("sources").size();reportProgress(0);
  Long64_t source=0,bin=0,responseBegin=0,responseEnd=0;double value=0,variance=0;std::string histogram,sourceJson;
  TTree contributions("source_contributions","Sparse response, measurement and run-dependent QA cells, including ROOT flow bins");
  contributions.Branch("source",&source);contributions.Branch("histogram",&histogram);
  contributions.Branch("bin",&bin);contributions.Branch("value",&value);contributions.Branch("variance",&variance);
  contributions.SetAutoFlush(-8*1024*1024);
  TTree sources("sources","Source identities, input hashes, producer metadata and normalization");
  sources.Branch("source",&source);sources.Branch("json",&sourceJson);
  sources.Branch("contribution_begin",&responseBegin);sources.Branch("contribution_end",&responseEnd);
  std::map<std::string,std::unique_ptr<TH1>> sums;
  std::map<std::string,std::unique_ptr<THnSparse>> sparse;
  Table events,photons,jets;Long64_t eventOffset=0;
  std::map<std::string,Table> triggerEvents,triggerPhotons,triggerJets;
  std::map<std::string,Long64_t> triggerOffsets;
  HistogramJob available=HistogramJob::object(),metadata=job.at("metadata");
  std::string grids;
  for(const auto& part:job.at("sources")) {
    sourceJson=part.at("metadata").dump();responseBegin=contributions.GetEntries();
    for(const auto& product:part.at("products")) {
      const std::string kind=product.at("kind");
      std::unique_ptr<TFile> in(TFile::Open(product.at("path").get<std::string>().c_str()));
      if(!in||in->IsZombie()) throw std::runtime_error("cannot open source product");
      if(metaText(*in,"profile")!=metadata.at("profile").get<std::string>()) throw std::runtime_error("product profile mismatch");
      if(metaText(*in,"selection")!=metadata.at("selection").dump()) throw std::runtime_error("product selection mismatch");
      if(kind=="response") {
        const auto text=metaText(*in,"grids"); if(!grids.empty()&&grids!=text) throw std::runtime_error("response grids mismatch");grids=text;
      }
      if(kind=="extensions") {
        const auto a=HistogramJob::parse(metaText(*in,"availability"));
        if(!available.empty()&&available!=a) throw std::runtime_error("mixed extension availability; use a uniform input release");available=a;
      }
      for(const auto& spec:job.at("inventory")) {
        reportProgress(source);
        if(spec.at("kind")!=kind) continue;
        const std::string key=spec.at("key"),path=spec.at("path");
        auto* obj=in->Get(key.c_str());
        if(!obj) {if(kind!="extensions") throw std::runtime_error("missing required histogram: "+key);continue;}
        if(auto* h=dynamic_cast<TH1*>(obj)) {
          if(h->GetSumw2N()!=h->GetNcells()) throw std::runtime_error("histogram has no explicit sumw2: "+key);
          for(int cell=0;cell<h->GetNcells();++cell) {
            const double v=h->GetBinContent(cell),e=h->GetBinError(cell);
            if(!std::isfinite(v)||!std::isfinite(e)) throw std::runtime_error("nonfinite histogram cell");
            const bool runQA=kind=="qa"&&(path=="qa/events/centrality_scaled_bits"||path=="qa/events/trigger_bits_raw"||path=="qa/events/centrality");
            if((kind=="response"||kind=="spectra"||runQA)&&(v!=0||e!=0)) {bin=cell;histogram=path;value=v;variance=h->GetSumw2()->At(cell);contributions.Fill();}
          }
          if(!sums.count(path)) {sums[path].reset(static_cast<TH1*>(h->Clone()));sums[path]->SetDirectory(nullptr);}
          else {sameAxes(*sums[path],*h);if(!sums[path]->Add(h)) throw std::runtime_error("histogram addition failed");}
        } else if(auto* h=dynamic_cast<THnSparse*>(obj)) {
          if(!sparse.count(path)) sparse[path].reset(static_cast<THnSparse*>(h->Clone()));else sparse[path]->Add(h);
        } else throw std::runtime_error("unsupported histogram object "+key);
      }
      if(kind=="spectra") {
        auto* ev=dynamic_cast<TTree*>(in->Get("table_events"));auto* ph=dynamic_cast<TTree*>(in->Get("table_photons"));auto* jt=dynamic_cast<TTree*>(in->Get("table_jets"));
        if(!ev||!ph||!jt) throw std::runtime_error("source event tables missing");
        events.append(*ev,"events",source,eventOffset,records);photons.append(*ph,"photons",source,eventOffset,records);jets.append(*jt,"jets",source,eventOffset,records);
        eventOffset+=ev->GetEntries();
        if(metadata.at("system")=="AuAu"&&metadata.at("family")=="data") {
          const auto selected=metadata.at("selection").at("systems").at("auau").at("data_trigger").get<std::string>();
          for(const auto& key:{selected}) {
            const std::string base=std::string("triggers/")+key+"/";
            auto* te=dynamic_cast<TTree*>(in->Get((base+"table_events").c_str()));
            auto* tp=dynamic_cast<TTree*>(in->Get((base+"table_photons").c_str()));
            auto* tj=dynamic_cast<TTree*>(in->Get((base+"table_jets").c_str()));
            if(!te||!tp||!tj) throw std::runtime_error("retained trigger event tables missing");
            auto* dest=directory(*out,std::string("records/triggers/")+key);
            triggerEvents[key].append(*te,"events",source,triggerOffsets[key],dest);
            triggerPhotons[key].append(*tp,"photons",source,triggerOffsets[key],dest);
            triggerJets[key].append(*tj,"jets",source,triggerOffsets[key],dest);
            triggerOffsets[key]+=te->GetEntries();
          }
        }
      }
    }
    responseEnd=contributions.GetEntries();sources.Fill();++source;reportProgress(source);
  }
  for(auto& item:sums) {const auto slash=item.first.rfind('/');directory(*out,item.first.substr(0,slash))->cd();item.second->Write(item.first.substr(slash+1).c_str());}
  for(auto& item:sparse) {const auto slash=item.first.rfind('/');directory(*out,item.first.substr(0,slash))->cd();item.second->Write(item.first.substr(slash+1).c_str());}
  records->cd();contributions.Write("",TObject::kOverwrite);sources.Write("",TObject::kOverwrite);
  events.write();photons.write();jets.write();
  for(auto& item:triggerEvents) item.second.write();
  for(auto& item:triggerPhotons) item.second.write();
  for(auto& item:triggerJets) item.second.write();
  out->cd();TObjString(metadata.at("profile").get<std::string>().c_str()).Write("profile");
  TObjString(metadata.at("selection").dump().c_str()).Write("selection");
  if(!grids.empty()) TObjString(grids.c_str()).Write("grids");
  metadata["histograms_written"]=sums.size()+sparse.size();metadata["availability"]=available;
  metadata["source_count"]=source;metadata["event_records"]=eventOffset;
  if(!triggerOffsets.empty()) {
    metadata["trigger_event_records"]=triggerOffsets;
    metadata["data_trigger"]=metadata.at("selection").at("systems").at("auau").at("data_trigger");
    metadata["trigger_retention"]="named_scaled_trigger_measurements_and_source_contributions";
  }
  directory(*out,"meta")->cd();TObjString(metadata.dump().c_str()).Write("package");
  TObjString(job.at("inventory").dump().c_str()).Write("inventory");
  // Stack-owned ROOT trees detach before closing the owning file.
  contributions.SetDirectory(nullptr);sources.SetDirectory(nullptr);out->Close();
}
}
