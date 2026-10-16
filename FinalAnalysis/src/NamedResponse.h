/* Selection of an already-produced response component. No energy operation,
   new normalization, event reduction or scientific prescription is applied. */
#pragma once
#include "FinalAnalysis.h"
#include "RenderMetadata.h"
#include "../../TreeToHists/src/CollectionContract.h"
#include <TFile.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TObjString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <limits>

namespace fa { namespace named_response {
using Json=nlohmann::json;
inline void require(bool ok,const std::string &why) {if(!ok)throw std::runtime_error("named response: "+why);}
inline bool pin(const std::string &text) {return text.size()==64 && text.find_first_not_of("0123456789abcdef")==std::string::npos;}
inline bool active(const SystemSettings &sys) {return !sys.responseComponent.empty();}
inline std::filesystem::path codeFile() {return (std::filesystem::path(__FILE__).parent_path()/"../../TreeToHists/src/CollectionContract.h").lexically_normal();}
inline std::string prefix(const SystemSettings &sys) {return active(sys)?"response_components/"+sys.responseComponent+"/":"";}
inline Json object(TFile &file,const char *key) {
  auto *value=dynamic_cast<TObjString*>(file.Get(key));require(value,std::string("missing ")+key);return Json::parse(value->GetString().Data());
}
/* Decode the existing persisted Variation format for compatibility checks.
   Values are never applied here: named resolution/scale draws are upstream. */
inline Json resolvedRecipe(const std::string &text,const std::string &system) {
  const auto items=photonjet_collection::variationItems(text);
  Json out={{"jet_scale",1.},{"jet_smear",0.},{"photon_scale",1.},{"photon_smear",0.},{"em_scale",1.},
    {"photon_resolution","off"},{"jet_resolution","off"},{"seed",0ULL}};
  bool only=false,off=false,ems=false;
  for(const auto &[key,value]:items) {
    if(key=="photon_selection") {require(value=="response_only" || value=="ems_selection_response","unknown selection fanout");only=value=="response_only";ems=value=="ems_selection_response";}
    else if(key=="centrality_weight") {require(value=="on" || value=="off","invalid centrality policy");off=value=="off";}
    else if(key=="photon_resolution") {
      require(value=="off" || value=="ppg18_nominal" || value=="ppg18_high" || value=="ppg18_low","unknown named photon resolution");out[key]=value;
    } else if(key=="jet_resolution") {
      require(value=="off" || value=="sphenix_pp2026_nominal" || value=="sphenix_pp2026_high" || value=="sphenix_pp2026_low","unknown named jet resolution");out[key]=value;
    } else if(key=="seed") {
      require(!value.empty() && value.find_first_not_of("0123456789")==std::string::npos,"invalid seed");out[key]=std::stoull(value);
    } else if(key=="jet_scale" || key=="photon_scale" || key=="em_scale" || key=="jet_smear" || key=="photon_smear") {
      size_t used=0;const double number=std::stod(value,&used);require(used==value.size() && std::isfinite(number) &&
        ((key=="jet_smear" || key=="photon_smear")?number>=0:number>0),"invalid scale/smear metadata");out[key]=number;
    } else require(false,"unknown recipe item: "+key);
  }
  require((ems || out.at("em_scale")==1.) && out.at("photon_scale")==1. && out.at("photon_smear")==0. && out.at("jet_smear")==0.,"unsupported component energy composition");
  if(ems) {
    photonjet_collection::requireEmsRecipe(text);
    require(out.at("photon_resolution")=="ppg18_nominal" && out.at("jet_resolution")=="sphenix_pp2026_nominal" &&
      out.at("jet_scale")==1. && !off,"EMS component requires nominal resolutions and weights");
    out["photon_selection"]="ems_selection_response";
  }
  if(only) {
    require(out.at("photon_resolution")!="off","response-only requires named photon resolution");out["photon_selection"]="response_only";
  }
  if(off) {
    const bool nominalResolution=only && out.at("photon_resolution")=="ppg18_nominal" &&
      out.at("jet_resolution")=="sphenix_pp2026_nominal";
    const bool pure=!only && out.at("photon_resolution")=="off" && out.at("jet_resolution")=="off";
    require(system=="AuAu" && out.at("jet_scale")==1. && (pure || nominalResolution),
      "centrality off requires AuAu and a pure weight alternative or nominal response-only resolutions");
    out["centrality_weight"]="off";
  }
  require(text.empty() || only || off || ems,"component is not nominal, response-only, or pure centrality off");return out;
}
inline void validateChoice(const SystemSettings &sys) {
  const bool any=active(sys) || !sys.responseRecipe.empty() || !sys.responseGridSha256.empty();if(!any)return;
  require(active(sys) && !sys.responseRecipe.empty() && pin(sys.responseGridSha256),"component, recipe and grid SHA256 are required together");
  require(sys.package=="PP" || sys.package=="AuAu","component system must be PP or AuAu");
  require(sys.responseComponent!="nominal" && sys.responseComponent.size()<=64 &&
    sys.responseComponent.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")==std::string::npos,"unsafe/nominal component name");
  resolvedRecipe(sys.responseRecipe,sys.package);
}
inline Json selection(const SystemSettings &sys) {
  validateChoice(sys);return {{"schema","PhotonJetNamedResponseSelectionV1"},{"component",sys.responseComponent},
    {"recipe",sys.responseRecipe},{"resolved_recipe",resolvedRecipe(sys.responseRecipe,sys.package)},
    {"grid_configuration_sha256",sys.responseGridSha256},{"scope","RESPONSE_COMPONENT_ONLY"}};
}
inline void sameAxis(const TAxis *a,const TAxis *b) {
  require(a && b && a->GetNbins()==b->GetNbins(),"component axis dimensions differ");
  for(int i=1;i<=a->GetNbins()+1;++i)require(a->GetBinLowEdge(i)==b->GetBinLowEdge(i),"component axis edges differ");
}
inline void declaredAxis(const TAxis *axis,const Json &edges) {
  require(axis && edges.is_array() && edges.size()==static_cast<size_t>(axis->GetNbins()+1),"declared response axis dimensions differ");
  for(int i=1;i<=axis->GetNbins()+1;++i)require(axis->GetBinLowEdge(i)==edges.at(i-1).get<double>(),"declared response axis edges differ");
}
/* Validate the persisted grid before reading its selected namespace. EMS uses
   the producer's pinned calibration-selection contract, never a resolution-only exemption. */
inline void verifyGridConfiguration(const Json &config,const Json &nominalSelection,const SystemSettings &sys) {
  const bool withEms=config.contains("ems_selection_contract");
  require((photonjet_collection::keys(config,{"schema","variants","limits"}) ||
    photonjet_collection::keys(config,{"schema","variants","limits","ems_selection_contract"})) && config.at("schema")=="PhotonJetNativeResponseGridSettingsV1" &&
    photonjet_render::sha256(config.dump())==sys.responseGridSha256,"grid configuration hash/schema differs");
  require(config.at("variants").is_array() && !config.at("variants").empty() && config.at("variants").size()<=16,"invalid component roster");
  require(photonjet_collection::keys(config.at("limits"),{"max_histogram_bytes","max_resident_bytes","max_output_bytes","max_sparse_bins_per_histogram","max_seconds"}),"component resource contract differs");
  for(const auto &[name,value]:config.at("limits").items())require(value.is_number_integer() && value.get<long long>()>0,"invalid component resource metadata");
  std::set<std::string> names,recipes;bool nominal=false,found=false,sawEms=false;
  for(const auto &variant:config.at("variants")) {
    require(photonjet_collection::keys(variant,{"name","variation"}) && variant.at("name").is_string() && variant.at("variation").is_string(),"invalid component declaration");
    const std::string name=variant.at("name"),recipe=variant.at("variation");
    require(!name.empty() && name.size()<=64 && name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")==std::string::npos && names.insert(name).second,"duplicate/unsafe component name");
    const auto resolved=resolvedRecipe(recipe,sys.package);
    if(photonjet_collection::emsRecipe(recipe)) {
      sawEms=true;require(withEms,"EMS grid is missing its selection contract");
      photonjet_collection::requireEmsContract(config.at("ems_selection_contract"));
      require(config.at("ems_selection_contract").at("selection")==nominalSelection,"EMS selection differs from package");
    }
    require(recipes.insert(resolved.dump()).second,"duplicate resolved component recipe");
    if(recipe.empty()){require(name=="nominal","nominal recipe has another label");nominal=true;}
    else require(name!="nominal","nominal component has a nonempty recipe");
    if(name==sys.responseComponent){require(resolved==resolvedRecipe(sys.responseRecipe,sys.package),"selected component recipe differs");found=true;}
  }
  require(nominal && found,"requested component is missing");
  require(sawEms==withEms,"EMS contract/variant roster differs");
}
/* The typed collection pins the package bytes; this checks the selected
   internal namespace against its declared grid, recipe, source and inventory. */
inline void verifyPackage(TFile &file,const SystemSettings &sys) {
  validateChoice(sys);require(active(sys),"missing explicit component");
  const auto meta=object(file,"meta/package");
  require(meta.at("schema")=="PhotonJetSampleV1" && meta.at("package_kind")=="FINAL" && meta.at("system")==sys.package && meta.at("family")=="photonjet","component needs a final signal-SIM package");
  const auto &histograms=meta.at("histograms");
  const auto grids=object(file,"grids");
  require(photonjet_collection::keys(grids,{"reco_pt","truth_pt","xj"}) && grids.at("reco_pt")==histograms.at("native_reco_response_edges_gev") &&
    grids.at("truth_pt")==histograms.at("native_truth_response_edges_gev") && grids.at("xj")==histograms.at("native_xj_edges"),"declared response grids differ from package settings");
  require(histograms.at("response_variation")=="" && !histograms.contains("measurement_replay") && !histograms.contains("ems_selection_contract") && !meta.contains("replay_provenance"),"a named component cannot overlay another response/replay/EMS recipe");
  const auto &normal=meta.at("normalization");require(normal.at("applied_in_histograms")==false,"component normalization was already applied");
  const auto &weights=meta.at("centrality_weights");
  require(weights.is_array() && weights.size()==(sys.package=="AuAu"?16U:0U),"component centrality map has the wrong system/dimensions");
  for(const auto &weight:weights)require(photonjet_collection::positive(weight),"component centrality map is not finite/positive");
  const auto &binding=meta.at("response_components");
  require(photonjet_collection::keys(binding,{"scope","configuration","ems_measurement_qa_fanout"}) && binding.at("scope")=="RESPONSE_COMPONENT_ONLY" && binding.at("ems_measurement_qa_fanout")=="UNQUALIFIED","invalid component scope");
  const auto &config=binding.at("configuration");
  verifyGridConfiguration(config,meta.at("selection"),sys);
  const auto &render=meta.at("render_selection");const auto &producer=render.at("producer_code");
  for(const char *name:{"TreeToHists.cc","FillPhotonJetResponse.C","EnergyVariation.h","CollectionContract.h","NativeResponseGrid.h"})
    require(producer.contains(name) && producer.at(name).is_string() && pin(producer.at(name)),"missing component producer pin");
  require(render.at("selection_sha256")==meta.at("selection").at("selection_sha256") && render.at("histogram_settings")==histograms,"component selector/axes metadata differs");
  auto *sources=dynamic_cast<TTree*>(file.Get("records/sources"));auto *cells=dynamic_cast<TTree*>(file.Get("records/source_contributions"));
  require(sources && cells && sources->GetEntries()==meta.at("source_count").get<Long64_t>() && sources->GetEntries()>0,"source-aware component records are missing");
  TTreeReader reader(sources);TTreeReaderValue<std::string> encoded(reader,"json");TTreeReaderValue<Long64_t> source(reader,"source"),begin(reader,"contribution_begin"),end(reader,"contribution_end");
  Long64_t expected=0,next=0;std::set<std::string> originals;
  while(reader.Next()) {
    const auto row=Json::parse(*encoded);const std::string identity=row.at("identity"),hash=row.at("sha256");
    require(*source==expected++ && *begin==next && *end>=*begin && *end<=cells->GetEntries() && pin(hash) && !identity.empty() && originals.insert(identity).second,"component source identity/contribution range differs");next=*end;
    const auto &origin=row.at("metadata").at("source");require(origin.value("logical_source_id",origin.value("input_file_sha256",std::string()))==identity,"component original-source identity differs");
    long long entry=0;for(const auto &range:row.at("entry_ranges")) {
      require(range.is_array() && range.size()==2 && range.at(0)==entry && range.at(1).get<long long>()>entry,"component original entry ranges differ");entry=range.at(1);
    }
    require(entry==row.at("reduced_entries") && entry==row.at("metadata").at("entries"),"component source coverage differs");
  }
  require(reader.GetEntryStatus()==TTreeReader::kEntryBeyondEnd && next==cells->GetEntries(),"component source records are incomplete");
  const auto inventory=object(file,"meta/inventory");std::map<std::string,Json> rows;
  for(const auto &row:inventory)require(rows.emplace(row.at("path"),row).second,"duplicate component inventory path");
  std::set<std::string> selected;const auto root=prefix(sys);
  for(const auto &[path,row]:rows)if(row.at("kind")=="response") {
    require(path.rfind("response/",0)==0 && path.find("..") == std::string::npos,"response inventory leaves its namespace");
    const std::string target=root+path;Json expected=row;expected["path"]=target;expected["kind"]="response__"+sys.responseComponent;
    require(rows.count(target) && rows.at(target)==expected,"component inventory is not the exact nominal response layout");selected.insert(target);
    auto *a=file.Get(path.c_str()),*b=file.Get(target.c_str());require(a && b,"component histogram missing");
    if(auto *h=dynamic_cast<TH1*>(a)) {
      auto *other=dynamic_cast<TH1*>(b);require(other && h->GetDimension()==other->GetDimension() && other->GetSumw2N()==other->GetNcells(),"component histogram type/Sumw2 differs");
      sameAxis(h->GetXaxis(),other->GetXaxis());if(h->GetDimension()>1)sameAxis(h->GetYaxis(),other->GetYaxis());if(h->GetDimension()>2)sameAxis(h->GetZaxis(),other->GetZaxis());
      const auto leaf=path.substr(path.rfind('/')+1);
      if(leaf=="xj_reco" || leaf=="xj_truth") {
        declaredAxis(h->GetXaxis(),grids.at(leaf=="xj_reco"?"reco_pt":"truth_pt"));declaredAxis(h->GetYaxis(),grids.at("xj"));
      }
    } else if(auto *h=dynamic_cast<THnSparse*>(a)) {
      auto *other=dynamic_cast<THnSparse*>(b);require(other && h->GetNdimensions()==other->GetNdimensions() && other->GetCalculateErrors(),"component sparse type/Sumw2 differs");
      for(int i=0;i<h->GetNdimensions();++i)sameAxis(h->GetAxis(i),other->GetAxis(i));
    } else require(false,"unsupported component object");
  }
  require(!selected.empty(),"empty response inventory");
  for(const auto &[path,row]:rows)if(path.rfind(root,0)==0)require(selected.count(path),"unrelated selected-component inventory entry");
}
inline void requirePinned(const Settings &s,const SystemSettings &sys,const std::string &sample,TFile &file) {
  if(!active(sys))return;
  const std::string key=sys.package+"/"+sample;
  require(!s.inputCollection.empty() && s.namedResponsePins.count(key),"component was not admitted through the typed collection");
  require(photonjet_render::sha256(file.GetName(),true)==s.namedResponsePins.at(key),"component package bytes changed after collection validation");
  verifyPackage(file,sys);
}
} }
