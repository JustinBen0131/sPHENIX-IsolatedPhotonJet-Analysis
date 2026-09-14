#ifndef RJ_SHOWER_FACTORIAL_V1_H
#define RJ_SHOWER_FACTORIAL_V1_H

#include "RJReplayFoundationV1.h"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace RJShowerFactorialV1
{
using namespace RJReplayFoundationV1;

enum class EnergySource : std::int32_t { CALIBRATED_TOWERINFO=0, RAWCLUSTER_MAP_VALUE=1 };
enum class Membership : std::int32_t { FULL_GRID=0, RAWCLUSTER_OWNED=1 };

struct Definition
{
  const char* name;
  EnergySource energy_source;
  Membership rectangular_membership;
  Membership moment_membership;
  double floor_gev;
};

inline const std::array<Definition,7>& definitions()
{
  static const std::array<Definition,7> values{{
      {"H70",EnergySource::CALIBRATED_TOWERINFO,Membership::FULL_GRID,Membership::RAWCLUSTER_OWNED,0.070},
      {"H0",EnergySource::CALIBRATED_TOWERINFO,Membership::FULL_GRID,Membership::RAWCLUSTER_OWNED,0.0},
      {"G70",EnergySource::CALIBRATED_TOWERINFO,Membership::FULL_GRID,Membership::FULL_GRID,0.070},
      {"G0",EnergySource::CALIBRATED_TOWERINFO,Membership::FULL_GRID,Membership::FULL_GRID,0.0},
      {"O70",EnergySource::CALIBRATED_TOWERINFO,Membership::RAWCLUSTER_OWNED,Membership::RAWCLUSTER_OWNED,0.070},
      {"O0",EnergySource::CALIBRATED_TOWERINFO,Membership::RAWCLUSTER_OWNED,Membership::RAWCLUSTER_OWNED,0.0},
      {"R70",EnergySource::RAWCLUSTER_MAP_VALUE,Membership::RAWCLUSTER_OWNED,Membership::RAWCLUSTER_OWNED,0.070}}};
  return values;
}

inline const Definition& definition(const std::string& name)
{
  for (const auto& value : definitions()) if (name == value.name) return value;
  throw std::invalid_argument("unknown shower definition: " + name);
}

inline std::string semanticText(const Definition& value)
{
  return std::string("RJ_SHOWER_DEFINITION_FACTORIAL_V1|") + value.name +
      "|energy=" + std::to_string(static_cast<int>(value.energy_source)) +
      "|sums=" + std::to_string(static_cast<int>(value.rectangular_membership)) +
      "|moments=" + std::to_string(static_cast<int>(value.moment_membership)) +
      "|floor_gev=" + (value.floor_gev > 0.0 ? "0.070000" : "0.000000") +
      "|grid=7x7|tower_quality=TowerInfo_get_isGood_only|center_excluded_from_cogx_numerator=1" +
      "|numeric=float32_PhotonClusterBuilder_row_major";
}

inline std::string semanticSha256(const Definition& value)
{
  static const std::array<std::string,7> cached=[]
  {
    std::array<std::string,7> hashes;
    const auto& values=definitions();
    for(std::size_t index=0;index<values.size();++index)
      hashes[index]=sha256Hex(semanticText(values[index]));
    return hashes;
  }();
  const auto& values=definitions();
  for(std::size_t index=0;index<values.size();++index)
    if(&value==&values[index])return cached[index];
  return sha256Hex(semanticText(value));
}

inline bool member(const ShowerCellRow& cell, Membership membership)
{
  return membership == Membership::FULL_GRID || cell.rawcluster_owned != 0;
}

inline bool selectedEnergy(const ShowerCellRow& cell,
                           const Definition& value,
                           float& energy)
{
  if (value.energy_source == EnergySource::CALIBRATED_TOWERINFO)
  {
    energy = static_cast<float>(cell.calibrated_energy);
    return cell.is_good != 0 && std::isfinite(energy) &&
        energy > static_cast<float>(value.floor_gev);
  }
  energy = static_cast<float>(cell.rawcluster_map_value);
  return cell.rawcluster_owned != 0 && cell.rawcluster_value_present != 0 &&
      std::isfinite(energy) && energy > static_cast<float>(value.floor_gev);
}

inline ShowerFeatureViewRow buildView(const Identity128& candidateId,
                                      const Definition& value,
                                      const std::vector<ShowerCellRow>& cells,
                                      double rawCenterEta,
                                      double rawCenterPhi,
                                      const std::array<double,4>& nativeEt)
{
  if (candidateId.isNull()) throw std::invalid_argument("null shower-view candidate identity");
  if (!std::isfinite(rawCenterEta) || !std::isfinite(rawCenterPhi))
    throw std::invalid_argument("nonfinite shower center");

  const int centerEta=static_cast<int>(std::floor(rawCenterEta));
  int centerPhi=static_cast<int>(std::floor(rawCenterPhi));
  while (centerPhi<0) centerPhi+=256;
  while (centerPhi>=256) centerPhi-=256;
  const float rawCenterEtaF=static_cast<float>(rawCenterEta);
  const float rawCenterPhiF=static_cast<float>(rawCenterPhi);
  const float cogEtaLocal=3.0F+
      (rawCenterEtaF-static_cast<float>(std::floor(rawCenterEtaF))-0.5F);
  const float cogPhiLocal=3.0F+
      (rawCenterPhiF-static_cast<float>(std::floor(rawCenterPhiF))-0.5F);
  auto deltaPhiIndex=[](int towerPhi,int referencePhi)
  {
    int delta=towerPhi-referencePhi;
    while(delta<-128)delta+=256;
    while(delta>127)delta-=256;
    return delta;
  };

  ShowerFeatureViewRow row;
  row.candidate_id=candidateId;
  row.definition_name=value.name;
  row.semantic_sha256=semanticSha256(value);
  row.definition_id=identityFromSha256(row.semantic_sha256);
  row.floor_gev=value.floor_gev;
  row.cog_eta=cogEtaLocal;
  row.cog_phi=cogPhiLocal;
  row.raw_center_eta=rawCenterEta;
  row.raw_center_phi=rawCenterPhi;
  row.center_eta_index=centerEta;
  row.center_phi_index=centerPhi;
  row.energy_source=static_cast<int>(value.energy_source);
  row.rectangular_membership=static_cast<int>(value.rectangular_membership);
  row.moment_membership=static_cast<int>(value.moment_membership);
  row.native_et1=nativeEt[0]; row.native_et2=nativeEt[1];
  row.native_et3=nativeEt[2]; row.native_et4=nativeEt[3];

  const int signPhi=cogPhiLocal>3.0?1:-1;
  float e11=0.0F,e33=0.0F,e32=0.0F,e35=0.0F;
  float momentEtaNumerator=0.0F,momentPhiNumerator=0.0F,momentDenominator=0.0F;
  float moment33EtaNumerator=0.0F,moment33PhiNumerator=0.0F,moment33Denominator=0.0F;
  for (const auto& cell : cells)
  {
    if (cell.candidate_id != candidateId) continue;
    const int i=cell.tower_eta_index-centerEta+3;
    const int j=deltaPhiIndex(cell.tower_phi_index,centerPhi)+3;
    if (i<0||i>6||j<0||j>6) continue;
    row.good_cell_count += cell.is_good != 0;
    row.owned_cell_count += cell.rawcluster_owned != 0;
    row.exact_zero_count += cell.is_good != 0 && cell.is_zero != 0;
    row.negative_count += cell.is_good != 0 && cell.is_negative != 0;
    row.nonfinite_count += cell.is_nonfinite != 0;

    float energy=std::numeric_limits<float>::quiet_NaN();
    if (!selectedEnergy(cell,value,energy)) continue;
    const bool sumMember=member(cell,value.rectangular_membership);
    const bool momentMember=member(cell,value.moment_membership);
    const int di=std::abs(i-3),dj=std::abs(j-3);
    if (sumMember)
    {
      ++row.active_sum_cell_count;
      if (i==3&&j==3) e11+=energy;
      if (di<=1&&dj<=1) e33+=energy;
      if (di<=1&&(j==3||j==3+signPhi)) e32+=energy;
      if (di<=1&&dj<=2) e35+=energy;
    }
    if (momentMember)
    {
      ++row.active_moment_cell_count;
      const float deta=static_cast<float>(i)-cogEtaLocal;
      const float dphi=static_cast<float>(j)-cogPhiLocal;
      momentDenominator+=energy;
      if (i!=3||j!=3)
      {
        momentEtaNumerator+=energy*deta*deta;
        momentPhiNumerator+=energy*dphi*dphi;
      }
      if (di<=1&&dj<=1)
      {
        moment33Denominator+=energy;
        if (i!=3||j!=3)
        {
          moment33EtaNumerator+=energy*deta*deta;
          moment33PhiNumerator+=energy*dphi*dphi;
        }
      }
    }
  }

  row.e11=e11; row.e33=e33; row.e32=e32; row.e35=e35;
  row.moment_eta_numerator=momentEtaNumerator;
  row.moment_phi_numerator=momentPhiNumerator;
  row.moment_denominator=momentDenominator;
  row.moment33_eta_numerator=moment33EtaNumerator;
  row.moment33_phi_numerator=moment33PhiNumerator;
  row.moment33_denominator=moment33Denominator;
  if (e33>0.0F) row.e11_over_e33=e11/e33;
  if (e35>0.0F) row.e32_over_e35=e32/e35;
  if (momentDenominator>0.0F)
  {
    row.weta_cogx=momentEtaNumerator/momentDenominator;
    row.wphi_cogx=momentPhiNumerator/momentDenominator;
  }
  if (moment33Denominator>0.0F)
  {
    row.weta33_cogx=moment33EtaNumerator/moment33Denominator;
    row.wphi33_cogx=moment33PhiNumerator/moment33Denominator;
  }
  row.finite_feature_state=
      std::isfinite(row.weta_cogx)&&std::isfinite(row.wphi_cogx)&&
      std::isfinite(row.weta33_cogx)&&std::isfinite(row.wphi33_cogx)&&
      std::isfinite(row.e11_over_e33)&&std::isfinite(row.e32_over_e35)&&
      std::all_of(nativeEt.begin(),nativeEt.end(),[](double item){return std::isfinite(item);});
  return row;
}

inline std::vector<float> modelFeatures(const ShowerFeatureViewRow& row,
                                        double photonEt,
                                        double vertexZ,
                                        double eta,
                                        double centrality,
                                        bool isAuAu)
{
  if (isAuAu)
  {
    // Canonical 14-feature Au+Au order.  Keep the 3x3 widths adjacent to
    // their cog widths and centrality last; this is the load-bearing model
    // input contract, not an internal convenience ordering.
    return {
        static_cast<float>(photonEt),static_cast<float>(row.weta_cogx),
        static_cast<float>(row.wphi_cogx),static_cast<float>(row.weta33_cogx),
        static_cast<float>(row.wphi33_cogx),static_cast<float>(vertexZ),
        static_cast<float>(eta),static_cast<float>(row.e11_over_e33),
        static_cast<float>(row.native_et1),static_cast<float>(row.native_et2),
        static_cast<float>(row.native_et3),static_cast<float>(row.native_et4),
        static_cast<float>(row.e32_over_e35),static_cast<float>(centrality)};
  }
  return {
      static_cast<float>(photonEt),static_cast<float>(row.weta_cogx),
      static_cast<float>(row.wphi_cogx),static_cast<float>(vertexZ),
      static_cast<float>(eta),static_cast<float>(row.e11_over_e33),
      static_cast<float>(row.native_et1),static_cast<float>(row.native_et2),
      static_cast<float>(row.native_et3),static_cast<float>(row.native_et4),
      static_cast<float>(row.e32_over_e35)};
}
} // namespace RJShowerFactorialV1

#endif
