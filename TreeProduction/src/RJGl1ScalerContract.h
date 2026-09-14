#ifndef RJGL1SCALERCONTRACT_H
#define RJGL1SCALERCONTRACT_H

#include <array>
#include <cstdint>
#include <limits>
#include <string>

namespace RJGl1ScalerContract
{
enum State : int { VALID=0, MISSING_PACKET=1, UNSUPPORTED_PACKET=2,
                   READ_ERROR=3, BAD_PACKET_STATUS=4 };
struct Snapshot
{
  std::array<std::uint64_t,64> raw{},live{},scaled{};
  std::uint64_t bco=std::numeric_limits<std::uint64_t>::max(), packet_status=0;
  std::int64_t physical_event=-1;
  int state=MISSING_PACKET,version=0;
};

// Gl1Packetv1 only overrides lValue(i,j); inherited getScaler returns zero.
// v2/v3 expose the counters through uint64_t getScaler. Preserve the high bit
// when converting the legacy signed API back to its stored uint64_t value.
template<class Packet> Snapshot read(const Packet* packet)
{
  Snapshot s;
  if (!packet) return s;
  const std::string type=packet->ClassName();
  s.version=type=="Gl1Packetv1" ? 1 : type=="Gl1Packetv2" ? 2 :
      type=="Gl1Packetv3" ? 3 : 0;
  if (!s.version) { s.state=UNSUPPORTED_PACKET; return s; }
  try
  {
    s.bco=packet->getBCO();
    s.physical_event=packet->getEvtSequence();
    s.packet_status=packet->getStatus();
    for (int i=0;i<64;++i)
    {
      auto counter=[&](int j)->std::uint64_t {
        return s.version==1 ? static_cast<std::uint64_t>(packet->lValue(i,j))
                            : packet->getScaler(i,j);
      };
      s.raw[i]=counter(0); s.live[i]=counter(1); s.scaled[i]=counter(2);
    }
    s.state=s.packet_status==0 ? VALID : BAD_PACKET_STATUS;
  }
  catch (...) { s.state=READ_ERROR; }
  return s;
}

inline bool identicalCounters(const Snapshot& a,const Snapshot& b)
{
  return a.state==b.state && a.version==b.version &&
    a.packet_status==b.packet_status && a.raw==b.raw &&
    a.live==b.live && a.scaled==b.scaled;
}

// These are evidence flags, never an automatic reset/wrap correction.
enum Discontinuity : unsigned { COUNTER_DECREASE=1, BCO_REGRESSION=2,
  EVENT_NONINCREASING=4, MISSING_BCO=8, MISSING_EVENT_ID=16 };
inline unsigned discontinuity(const Snapshot* previous,const Snapshot& next)
{
  unsigned flags=0;
  constexpr auto unknown=std::numeric_limits<std::uint64_t>::max();
  if (next.bco==unknown) flags|=MISSING_BCO;
  if (next.physical_event<0) flags|=MISSING_EVENT_ID;
  if (!previous) return flags;
  if (previous->bco!=unknown && next.bco!=unknown && next.bco<previous->bco)
    flags|=BCO_REGRESSION;
  if (previous->physical_event>=0 && next.physical_event>=0 &&
      next.physical_event<=previous->physical_event) flags|=EVENT_NONINCREASING;
  if (previous->state==VALID && next.state==VALID)
    for (int i=0;i<64;++i)
      if (next.raw[i]<previous->raw[i] || next.live[i]<previous->live[i] ||
          next.scaled[i]<previous->scaled[i]) flags|=COUNTER_DECREASE;
  return flags;
}
}
#endif
