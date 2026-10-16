#pragma once
/* Producer provenance, separate from effective scientific definitions. The
   current binary never replaces a missing input binary's ancestry. */
#include "ScientificTrace.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TObjString.h>
#include <memory>

namespace photonjet_science {
inline Json compiledProducerWitness() {
#ifdef PHOTONJET_NATIVE_SOURCE_WITNESS
  const auto value=Json::parse(PHOTONJET_NATIVE_SOURCE_WITNESS);
  require(value.value("schema","")=="NativeProducerSourceWitnessV1" && value.value("state","")=="recorded" &&
    value.contains("source_inventory") && !value.at("source_inventory").empty(),"invalid compiled producer witness");
  return value;
#else
  return {{"schema","NativeProducerSourceWitnessV1"},{"state","unknown"},{"reason","binary has no compiled full source witness"}};
#endif
}
inline Json producerRecord(TDirectory& input,const char* key,const Json& missing) {
  TObject* saved=input.Get(key);if(!saved)saved=input.Get((std::string("meta/")+key).c_str());
  if(!saved)return missing;
  auto* text=dynamic_cast<TObjString*>(saved);require(text,"producer provenance must be a native JSON string");
  return Json::parse(text->GetString().Data());
}
inline Json missingProducerWitness() {
  return {{"schema","NativeProducerSourceWitnessV1"},{"state","unknown"},
    {"reason","input has no compiled source witness"}};
}
inline Json missingProducerAncestry() {
  return {{"schema","NativeProducerAncestryV1"},{"state","unknown"},
    {"reason","input ancestry was not recorded"},{"inputs",Json::array()}};
}
inline Json consolidatedProducerAncestry(const Json& consolidation) {
  require(consolidation.is_object() && consolidation.contains("chunks") &&
    consolidation.at("chunks").is_array(),"consolidation provenance lacks its exact chunk roster");
  Json inputs=Json::array();
  for(const auto& chunk:consolidation.at("chunks")) {
    const auto& metadata=chunk.at("original_meta");
    require(metadata.is_object(),"consolidated original metadata must be an object");
    inputs.push_back({{"source_path",chunk.at("path")},{"source_sha256",chunk.at("sha256")},
      {"entry_range",{chunk.at("entry_begin"),chunk.at("entry_end")}},
      {"witness",metadata.value("native_source_witness",missingProducerWitness())},
      {"ancestry",metadata.value("native_source_ancestry",missingProducerAncestry())}});
  }
  return {{"schema","NativeProducerAncestryV1"},{"state",inputs.empty()?"unknown":"recorded"},
    {"operation","eventtree.consolidate/v1"},{"inputs",inputs},
    {"qualification","REQUIRES_RETAINED_EXECUTABLE_RELEASE_AND_EVERY_PARENT"}};
}
inline Json producerInput(TDirectory& input) {
  auto witness=producerRecord(input,"native_source_witness",missingProducerWitness());
  auto ancestry=producerRecord(input,"native_source_ancestry",missingProducerAncestry());
  const auto consolidated=producerRecord(input,"consolidation",Json());
  if(!consolidated.is_null()) {
    const auto actual=consolidatedProducerAncestry(consolidated);
    if(ancestry.value("operation",std::string())=="eventtree.consolidate/v1")
      require(ancestry==actual,"consolidation ancestry differs from its retained original chunk roster");
    else {
      // Older consolidation copied the first chunk's metadata header. That
      // header cannot stand in for either the consolidator or other chunks.
      witness=missingProducerWitness();
      witness["reason"]="legacy consolidator did not capture its own compiled source witness";
    }
    ancestry=actual;
  }
  return {{"source_path",input.GetFile()?input.GetFile()->GetName():input.GetName()},
    {"witness",witness},{"ancestry",ancestry}};
}
inline Json producerInputs(const std::vector<std::string>& paths) {
  Json values=Json::array();
  for(const auto& path:paths) {
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(),"READ"));require(file && !file->IsZombie(),"cannot read producer ancestry: "+path);
    values.push_back(producerInput(*file));
  }
  return values;
}
inline void writeProducerEvidence(TDirectory& output,const Json& inputs) {
  require(inputs.is_array(),"producer ancestry requires the exact input roster");
  TDirectory::TContext context(&output);
  TObjString(compiledProducerWitness().dump().c_str()).Write("native_source_witness",TObject::kOverwrite);
  const Json ancestry={{"schema","NativeProducerAncestryV1"},{"state",inputs.empty()?"unknown":"recorded"},
    {"inputs",inputs},{"qualification","REQUIRES_RETAINED_EXECUTABLE_RELEASE_AND_EVERY_PARENT"}};
  TObjString(ancestry.dump().c_str()).Write("native_source_ancestry",TObject::kOverwrite);
}
}
