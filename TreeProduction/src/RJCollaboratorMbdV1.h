#ifndef RJ_COLLABORATOR_MBD_V1_H
#define RJ_COLLABORATOR_MBD_V1_H

#include <cmath>
#include <limits>

// Capture the calibrated MbdPmtContainer values, without changing selection.
// Row has value-owned vectors. Call once on a fresh EventBundle before any
// reco-dependent serialization. Validity concerns charge and channel identity,
// not timing, fired-PMT selection, or calibration-quality certification.
namespace RJCollaboratorMbdV1
{
template<class Row, class Container, class Geometry>
void capture(Row& row, Container* pmts, Geometry* geometry)
{
  row.mbd_pmt_id.clear();
  row.mbd_pmt_arm.clear();
  row.mbd_pmt_charge.clear();
  row.mbd_pmt_charge_valid.clear();
  row.mbd_pmt_available = pmts ? 1 : 0;
  if (!pmts) return;
  const int count = pmts->get_npmt();
  if (count < 0 || count > 128)
  {
    row.mbd_pmt_available = -1;
    return;
  }
  for (int channel = 0; channel < count; ++channel)
  {
    const auto* pmt = pmts->get_pmt(channel);
    const double charge = pmt ? pmt->get_q()
        : std::numeric_limits<double>::quiet_NaN();
    const int arm = geometry ? geometry->get_arm(channel) : -1;
    row.mbd_pmt_id.push_back(channel); // container channel; never compact holes
    row.mbd_pmt_arm.push_back(arm == 0 || arm == 1 ? arm : -1);
    row.mbd_pmt_charge.push_back(charge); // keep zero and negative values
    row.mbd_pmt_charge_valid.push_back(
        pmt && pmt->get_pmt() == channel && std::isfinite(charge) ? 1 : 0);
  }
}
}
#endif
