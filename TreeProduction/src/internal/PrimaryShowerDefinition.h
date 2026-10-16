/* Bind primary shower observables to arguments captured at the builder call.
   An old builder without these witnesses stays unknown. No label or model
   filename is used to supply a missing floor or membership choice. */
#ifndef PHOTONJET_PRIMARY_SHOWER_DEFINITION_H
#define PHOTONJET_PRIMARY_SHOWER_DEFINITION_H
#include "Types.h"
#include "../../../TreeToHists/src/ScientificTrace.h"
namespace photonjet_science {
inline Json primaryShowerDefinition(const photonjet::PhotonRecord& p) {
  if(!std::isfinite(p.primaryShowerFloorGeV) || p.primaryShowerFloorGeV<0 ||
     p.primaryShowerEnergy<0 || p.primaryShowerEnergy>1 || p.primaryShowerSums<0 || p.primaryShowerSums>1 ||
     p.primaryShowerMoments<0 || p.primaryShowerMoments>1 || p.primaryShowerAlgorithm!=1)
    return unknown("primary builder did not capture a recognized complete calculation witness");
  return {{"schema","PrimaryShowerCalculationV1"},{"algorithm","PhotonClusterBuilder.calculate_shower_shapes/v1"},
      {"tower_floor_gev",real(p.primaryShowerFloorGeV)},{"floor_comparison","energy > floor"},
      {"energy_source",p.primaryShowerEnergy==0?"calibrated":"cluster_map"},
      {"rectangular_membership",p.primaryShowerSums==0?"grid":"owned"},
      {"moment_membership",p.primaryShowerMoments==0?"grid":"owned"},
      {"floor_precision","binary32"}};
}
inline Json primaryBranchDefinitions(const std::set<std::string>& observed) {
  Json calculations=Json::array();for(const auto& value:observed)calculations.push_back(Json::parse(value));
  if(calculations.empty())calculations.push_back(unknown("no primary calculation witness was observed"));
  Json definitions=Json::object();
  for(const auto* name:photonjet::PhotonRecord::primaryShapeNames) {
    if(std::string(name)=="bdt_score")continue; // model function has its own definition
    definitions[name]={{"calculation",calculations},{"observable",name}};
  }
  definitions["e11_over_e33"]={{"calculation",calculations},{"observable","e11/e33"},{"division","binary64; denominator > 0"}};
  definitions["e32_over_e35"]={{"calculation",calculations},{"observable","e32/e35"},{"division","binary64; denominator > 0"}};
  return definitions;
}
inline Json namedShowerDefinition(const photonjet::ShowerShapeRecord& p) {
  photonjet::PhotonRecord witness;
  witness.primaryShowerFloorGeV=p.calculationFloorGeV;
  witness.primaryShowerEnergy=p.energySource;witness.primaryShowerSums=p.rectangularMembership;
  witness.primaryShowerMoments=p.momentMembership;witness.primaryShowerAlgorithm=p.calculationAlgorithm;
  return primaryShowerDefinition(witness);
}
inline Json namedShowerDefinitions(const std::map<std::string,std::set<std::string>>& observed) {
  Json result=Json::object();
  for(const auto& row:observed) {
    Json values=Json::array();for(const auto& value:row.second)values.push_back(Json::parse(value));
    result[row.first]={{"calculation",values},{"model_input_precision","binary32"}};
  }
  return result;
}
}  // namespace photonjet_science
#endif
