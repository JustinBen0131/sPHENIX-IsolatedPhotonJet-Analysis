#ifndef RJ_GL1_TRIGGER_CONTRACT_H
#define RJ_GL1_TRIGGER_CONTRACT_H

#include <cstdint>
#include <string>

// Direct packet decisions, independent of event retention and of trigger-name
// helpers. Never interpret a missing/unsupported accessor's zero as observed.
namespace RJGl1TriggerContract
{
constexpr int kVersion = 1;
enum State : int { VALID=0, MISSING_PACKET=1, UNSUPPORTED_PACKET=2,
                   READ_ERROR=3, BAD_PACKET_STATUS=4, NOT_DATA=5 };
enum Availability : unsigned { LIVE=1, SCALED=2, LEGACY=4 };
struct Decisions
{
  std::uint64_t legacy=0, live=0, scaled=0, packet_status=0;
  unsigned available=0;
  int version=0, state=MISSING_PACKET;
};

template<class Packet> Decisions read(const Packet* packet, bool isSimulation=false)
{
  Decisions d;
  if (isSimulation) { d.state=NOT_DATA; return d; }
  if (!packet) return d;
  const std::string type=packet->ClassName();
  d.version=type=="Gl1Packetv1" ? 1 : type=="Gl1Packetv2" ? 2 :
            type=="Gl1Packetv3" ? 3 : 0;
  if (!d.version) { d.state=UNSUPPORTED_PACKET; return d; }
  try
  {
    d.packet_status=packet->getStatus();
    d.legacy=packet->getTriggerVector(); d.available=LEGACY;
    // v1 does not override live/scaled accessors. Preserve its legacy word
    // without claiming either modern decision word was measured.
    if (d.version>=2)
    {
      d.live=packet->getLiveVector(); d.scaled=packet->getScaledVector();
      d.available|=LIVE|SCALED;
    }
    d.state=d.packet_status ? BAD_PACKET_STATUS :
            (d.version>=2 ? VALID : UNSUPPORTED_PACKET);
  }
  catch (...) { d=Decisions{}; d.state=READ_ERROR; }
  return d;
}

inline bool scaledFired(const Decisions& d,int bit)
{
  return d.state==VALID && (d.available&SCALED) && bit>=0 && bit<64 &&
         ((d.scaled>>bit)&std::uint64_t{1});
}

// Name lookup is only for existing per-trigger histograms. Tree capture never
// depends on it. Enumerating names avoids the base API's ambiguous bit-0
// fallback for an unknown name. All 64 output bits are supported.
template<class Packet,class RunInfo>
bool firedByName(const Packet* packet,const RunInfo* info,const std::string& name)
{
  if (!info || name.empty() || name=="unknown" ||
      std::string(info->ClassName())!="TriggerRunInfov1") return false;
  int bit=-1;
  for(int i=0;i<64;++i)
    if(info->getTriggerName(i)==name)
    {
      if(bit>=0) return false; // ambiguous name, not a guessed decision
      bit=i;
    }
  return scaledFired(read(packet),bit);
}

template<class Event> void capture(Event& event,const Decisions& d)
{
  event.trigger_capture_version=kVersion;
  event.trigger_bits=d.legacy;
  event.live_trigger_bits=d.live;
  event.scaled_trigger_bits=d.scaled;
  event.trigger_decision_state=d.state;
  event.trigger_decision_available=d.available;
  event.trigger_packet_version=d.version;
  event.trigger_packet_status=d.packet_status;
}
}
#endif
