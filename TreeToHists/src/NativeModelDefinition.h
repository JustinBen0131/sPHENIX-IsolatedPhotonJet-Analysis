/* Read the scientific function of a retained ROOT RBDT, without evaluating,
   fitting or rewriting it. The file digest belongs in provenance; only the
   model's numeric trees, base response, transform and feature ordering enter
   scientific identity. Unknown model representations are refused. */
#ifndef PHOTONJET_NATIVE_MODEL_DEFINITION_H
#define PHOTONJET_NATIVE_MODEL_DEFINITION_H
#include "ScientificTrace.h"
#include <TBufferJSON.h>
#include <TFile.h>
#include <TMVA/RBDT.hxx>
#include <memory>
#include <mutex>

namespace photonjet_science {
inline Json rbdtDefinition(const Json& saved,const Json& featureOrder) {
  require(featureOrder.is_array() && !featureOrder.empty(),"RBDT needs its actual ordered features");
  std::set<std::string> features;
  for(const auto& feature:featureOrder)require(feature.is_string() && !feature.get<std::string>().empty() &&
      features.insert(feature.get<std::string>()).second,"RBDT feature order must be named and unique");
  const std::set<std::string> integers={"fRootIndices","fCutIndices","fLeftIndices","fRightIndices","fTreeNumbers"};
  const std::set<std::string> reals={"fCutValues","fResponses","fBaseResponses"};
  Json function={{"schema","ROOTRBDTScientificFunctionV1"},{"input_precision","binary32"},
                 {"feature_order",featureOrder},{"arithmetic","ROOT::TMVA::Experimental::RBDT_v1"}};
  for(const auto& name:integers) {
    require(saved.contains(name) && saved.at(name).is_array(),"RBDT integer array absent: "+name);
    for(const auto& v:saved.at(name))require(v.is_number_integer(),"RBDT index is not integral");
    function[name]=saved.at(name);
  }
  for(const auto& name:reals) {
    require(saved.contains(name) && saved.at(name).is_array(),"RBDT value array absent: "+name);
    function[name]=Json::array();
    for(const auto& value:saved.at(name))function[name].push_back(real(static_cast<double>(value.get<float>())));
  }
  require(saved.contains("fBaseScore") && saved.at("fLogistic").is_boolean(),"RBDT response transform absent");
  function["fBaseScore"]=real(static_cast<double>(saved.at("fBaseScore").get<float>()));
  function["fLogistic"]=saved.at("fLogistic");
  const auto n=function.at("fCutIndices").size();
  require(function.at("fCutValues").size()==n && function.at("fLeftIndices").size()==n && function.at("fRightIndices").size()==n,
          "RBDT node arrays disagree");
  for(const auto& value:function.at("fCutIndices"))require(value.get<std::size_t>()<featureOrder.size(),"RBDT feature slot outside contract");
  const auto& roots=saved.at("fRootIndices");
  require(saved.at("fTreeNumbers").size()==roots.size(),"RBDT tree/output roster differs");
  require(!saved.at("fBaseResponses").empty(),"RBDT has no output channel");
  std::vector<int> state(n,0);
  std::function<void(int)> visit=[&](int node) {
    require(node>=0 && std::size_t(node)<n,"RBDT root/child node outside array");
    require(state[node]!=1,"RBDT decision graph is cyclic");if(state[node]==2)return;
    state[node]=1;
    for(const auto* side:{"fLeftIndices","fRightIndices"}) {
      const int child=saved.at(side).at(node).get<int>();
      if(child>0)visit(child);
      else require(-static_cast<long long>(child)<static_cast<long long>(saved.at("fResponses").size()),"RBDT leaf outside response array");
    }
    state[node]=2;
  };
  for(std::size_t tree=0;tree<roots.size();++tree) {
    const int channel=saved.at("fTreeNumbers").at(tree).get<int>();
    require(channel>=0 && std::size_t(channel)<saved.at("fBaseResponses").size(),"RBDT output channel outside contract");
    visit(roots.at(tree).get<int>());
  }
  for(const auto value:state)require(value==2,"RBDT retains an unreachable numerical node");
  return function;
}
inline Json loadedRBDTDefinition(const TMVA::Experimental::RBDT& model,const Json& featureOrder) {
  // 9 decimal digits round-trip every binary32 value used by ROOT RBDT.
  // TBufferJSON's defaults are restored to avoid altering another writer.
  // ROOT keeps the supplied format pointer. Restoring a temporary string's
  // c_str() would leave global serialization dangling after this function.
  static std::mutex mutex;std::lock_guard<std::mutex> lock(mutex);
  // Keep these few interned strings alive through ROOT shutdown or unloading
  // this helper library; ROOT's global format pointers outlive that library.
  static auto* formats=new std::set<std::string>;
  const char* oldFloat=formats->emplace(TBufferJSON::GetFloatFormat()).first->c_str();
  const char* oldDouble=formats->emplace(TBufferJSON::GetDoubleFormat()).first->c_str();
  struct Restore {const char *f,*d;~Restore(){TBufferJSON::SetFloatFormat(f);TBufferJSON::SetDoubleFormat(d);}} restore{oldFloat,oldDouble};
  TBufferJSON::SetFloatFormat("%.9g");TBufferJSON::SetDoubleFormat("%.17g");
  return rbdtDefinition(Json::parse(TBufferJSON::ToJSON(&model).Data()),featureOrder);
}
inline Json readRBDTDefinition(const std::string& path,const std::string& key,const Json& featureOrder) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(),"READ"));
  require(file && !file->IsZombie(),"retained RBDT cannot be opened");
  auto* model=file->Get<TMVA::Experimental::RBDT>(key.c_str());
  require(model,"retained RBDT key is absent");
  return loadedRBDTDefinition(*model,featureOrder);
}
}  // namespace photonjet_science
#endif
