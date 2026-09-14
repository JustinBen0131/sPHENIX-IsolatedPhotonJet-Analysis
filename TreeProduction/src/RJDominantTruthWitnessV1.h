#ifndef RJ_DOMINANT_TRUTH_WITNESS_V1_H
#define RJ_DOMINANT_TRUTH_WITNESS_V1_H

#include <cmath>
#include <cstdint>
#include <limits>

// Capture the accepted energy-dominant-primary association, not a nearest-DR
// substitute and not a derived training label. The surrounding candidate row
// supplies the source/event/cluster identity. No ID, isolation or pT cut lives
// here. A missing evaluator is never an ordinary background label.
namespace RJDominantTruthWitnessV1
{
enum State : std::int32_t
{
  NOT_SIMULATION = -1, EVALUATOR_UNAVAILABLE = 0, NO_PRIMARY = 1,
  VALID_PRIMARY = 2, INVALID_PRIMARY = 3
};
enum EvaluatorMode : std::int32_t
{
  NOT_APPLICABLE = -1, UNAVAILABLE = 0, TOWERINFO = 1, LEGACY_RAWTOWER = 2
};
struct Witness
{
  std::int32_t state=NOT_SIMULATION, evaluator_mode=NOT_APPLICABLE;
  std::int32_t track_id=-1, pid=0, barcode=-1, embedding_id=0;
  double energy_contribution=std::numeric_limits<double>::quiet_NaN();
  std::int32_t vertex_id=-1;
};

template<class Evaluator, class Cluster, class TruthInfo>
Witness capture(bool simulation, Evaluator* evaluator, Cluster* cluster,
                TruthInfo* truth, EvaluatorMode mode)
{
  Witness result;
  if (!simulation) return result;
  result.state=EVALUATOR_UNAVAILABLE;
  result.evaluator_mode=mode;
  if (!evaluator || !cluster || !truth || mode==UNAVAILABLE ||
      mode==NOT_APPLICABLE || !evaluator->has_reduced_node_pointers())
    return result;
  auto* primary=evaluator->max_truth_primary_particle_by_energy(cluster);
  result.state=NO_PRIMARY;
  if (!primary) return result;
  result.track_id=primary->get_track_id();
  result.pid=primary->get_pid();
  result.barcode=primary->get_barcode();
  result.vertex_id=primary->get_vtx_id();
  result.energy_contribution=evaluator->get_energy_contribution(cluster,primary);
  result.state=INVALID_PRIMARY;
  if (result.track_id<=0 || result.pid==0 ||
      !std::isfinite(result.energy_contribution) || result.energy_contribution<0.)
    return result;
  result.embedding_id=truth->isEmbeded(result.track_id);
  result.state=VALID_PRIMARY;
  return result;
}
}
#endif
