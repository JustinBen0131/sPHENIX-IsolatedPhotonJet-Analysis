/* Retain exact producer contexts through candidate-cache source remapping.
   Contexts are measured definitions, never inferred from a source digest. */
#pragma once
#include "ScientificTrace.h"
#include <TTree.h>
#include <TNamed.h>
#include <TList.h>
namespace photonjet_science {
inline Json candidateContexts(const TTree& tree) {
  auto* saved=dynamic_cast<TNamed*>(const_cast<TTree&>(tree).GetUserInfo()->FindObject("physics_candidate_contexts"));
  return saved?Json::parse(saved->GetTitle()):Json::object();
}
inline void candidateContexts(TTree& tree,const Json& rows) {
  require(rows.is_object(),"candidate scientific source inventory is not an object");
  auto* list=tree.GetUserInfo();if(auto* old=list->FindObject("physics_candidate_contexts"))list->Remove(old);
  list->Add(new TNamed("physics_candidate_contexts",rows.dump().c_str()));
}
inline Json mergeCandidateContexts(const Json& old,const Json& incoming,long long source,
                                   const std::vector<Long64_t>* mapping) {
  Json out=old;
  for(const auto& [id,context]:incoming.items()) {
    const auto slot=std::stoll(id);require(slot>=0,"negative context source slot");
    require(!mapping||static_cast<std::size_t>(slot)<mapping->size(),"context outside source map");
    const std::string mapped=std::to_string(mapping?mapping->at(slot):source);
    require(!out.contains(mapped)||out.at(mapped)==context,"source context changed during merge");out[mapped]=context;
  }
  return out;
}
}
