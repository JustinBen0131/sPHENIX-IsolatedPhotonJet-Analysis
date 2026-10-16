#pragma once
#include "CandidateScientificContext.h"
#include "MeasurementScientificDefinitions.h"
namespace photonjet_science {
inline std::vector<Json> replayScientificContexts(TTree& events,const photonjet::Profile& profile,
    long long source,const Json& variation,const Json& axes,double normalization,bool timing,
    const Json* centralityWeights=nullptr,bool centralityOff=false) {
  const auto saved=candidateContexts(events);std::vector<Json> contexts;
  const auto include=[&](Json c) {
    if(!c.is_object()||c.value("profile_witness",Json())!=measurementProfileWitness(profile)) {
      // The numerical replay may still be allowed by its independent source
      // validator; no unsupported ancestry becomes a complete receipt.
      c=Json{{"system",profile.system==photonjet::System::pp?"PP":"AuAu"},{"simulation",profile.simulation},
        {"source_channel",{{"sample","unavailable"}}}};
    }
    c["energy_variation"]=variation;c["measurement_axes"]=axes;
    c["replay_normalization"]=normalization;
    if(timing)c["timing_selection"]={{"rule","finite(time_ns) and abs(time_ns-center_ns)<half_width_ns"},
      {"center_ns",-2.},{"half_width_ns",6.},{"time_definition",child(c,"/cluster_time_reference")}};
    if(centralityWeights)c["centrality_weights"]=*centralityWeights;
    if(centralityOff && c.contains("centrality_weights"))for(auto& value:c["centrality_weights"])value=1.;
    contexts.push_back(c);
  };
  if(source>=0)include(saved.value(std::to_string(source),Json()));
  else for(const auto& c:saved)include(c);
  if(contexts.empty())include(Json());
  return contexts;
}
}
