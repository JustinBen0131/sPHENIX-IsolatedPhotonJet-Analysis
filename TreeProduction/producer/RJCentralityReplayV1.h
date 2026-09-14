#ifndef RJ_CENTRALITY_REPLAY_V1_H
#define RJ_CENTRALITY_REPLAY_V1_H

#include <cmath>
#include <limits>

// Capture the explicitly bound centrality nodes, independently of the legacy
// collaborator charge-sum node. In embedding these need not be the same node.
// Calibrated q is in MBD get_q units; time is ns; MbdOut z is cm.
namespace RJCentralityReplayV1 {
template<class Row, class Pmts, class Out, class MinimumBias>
void capture(Row& row, Pmts* pmts, Out* out, MinimumBias* mb) {
  row.centrality_replay_version=1;
  row.centrality_pmt_available=pmts ? 1 : 0;
  row.centrality_pmt_id.clear(); row.centrality_pmt_charge.clear();
  row.centrality_pmt_time.clear(); row.centrality_pmt_valid.clear();
  row.centrality_pmt_selected.clear();
  row.centrality_mbd_z=std::numeric_limits<double>::quiet_NaN();
  row.centrality_mbd_z_valid=0;
  row.centrality_mbd_event=-1; row.centrality_mbd_clock=-1;
  row.centrality_mbd_femclock=-1;
  row.centrality_mb_decision=mb ? (mb->isAuAuMinimumBias() ? 1 : 0) : -1;
  row.centrality_selected_charge=std::numeric_limits<double>::quiet_NaN();
  row.centrality_inputs_valid=0;
  if(out) {
    row.centrality_mbd_z=out->get_zvtx();
    row.centrality_mbd_z_valid=std::isfinite(row.centrality_mbd_z) ? 1 : 0;
    row.centrality_mbd_event=out->get_evt();
    row.centrality_mbd_clock=out->get_clock();
    row.centrality_mbd_femclock=out->get_femclock();
  }
  if(!pmts) return;
  if(pmts->get_npmt()!=128) { row.centrality_pmt_available=-1; return; }
  bool valid=true;
  float selected_charge=0;
  for(int i=0;i<128;++i) {
    const auto* pmt=pmts->get_pmt(i);
    const float q=pmt ? pmt->get_q() : std::numeric_limits<float>::quiet_NaN();
    const float t=pmt ? pmt->get_time() : std::numeric_limits<float>::quiet_NaN();
    const bool good=pmt && pmt->get_pmt()==i && std::isfinite(q) && std::isfinite(t);
    const bool selected=good && q>=0.5f && std::fabs(t)<=25.0f;
    row.centrality_pmt_id.push_back(i);
    row.centrality_pmt_charge.push_back(q); row.centrality_pmt_time.push_back(t);
    row.centrality_pmt_valid.push_back(good ? 1 : 0);
    row.centrality_pmt_selected.push_back(selected ? 1 : 0);
    valid=valid && good;
    if(selected) selected_charge+=q;
  }
  if(valid) row.centrality_selected_charge=selected_charge;
  row.centrality_inputs_valid=(valid && row.centrality_mbd_z_valid && mb) ? 1 : 0;
}
} // namespace RJCentralityReplayV1
#endif
