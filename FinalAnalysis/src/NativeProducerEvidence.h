#pragma once
/* Preserve the exact inputs consumed by this invocation. Provenance is not
   scientific identity; a new executable cannot qualify an older input. */
#include "../../TreeToHists/src/ScientificProducerEvidence.h"
namespace fa {
inline nlohmann::json& nativeInputs() {
  static thread_local nlohmann::json inputs=nlohmann::json::array();
  return inputs;
}
inline void resetNativeInputs() { nativeInputs()=nlohmann::json::array(); }
inline void observeNativeInput(TDirectory& input) {
  const auto record=photonjet_science::producerInput(input);
  for(const auto& previous:nativeInputs()) if(previous==record) return;
  nativeInputs().push_back(record);
}
inline nlohmann::json nativeProducerWitness() { return photonjet_science::compiledProducerWitness(); }
inline nlohmann::json nativeProducerAncestry() {
  return {{"schema","NativeProducerAncestryV1"},{"state",nativeInputs().empty()?"unknown":"recorded"},
          {"inputs",nativeInputs()}};
}
inline void writeNativeProducerWitness(TDirectory& output) {
  photonjet_science::writeProducerEvidence(output,nativeInputs());
}
}
