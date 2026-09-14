#ifndef RJTRIGGERSCALERSV1_H
#define RJTRIGGERSCALERSV1_H

#include "RJGl1ScalerContract.h"
#include <TDirectory.h>
#include <TNamed.h>
#include <TTree.h>
#include <array>
#include <string>

namespace RJTriggerScalersV1
{
// Lossless counter snapshots, run-length encoded only when all 192 counters,
// validity and packet implementation agree. Each observed event receives its
// source-local snapshot_id; IDs from different physical sources never alias.
// Cumulative counters are NOT luminosities and must never be summed over rows.
class Writer
{
 public:
  bool initialize(TDirectory* directory,int run,int segment,
                  const std::string& inputHash,const std::string& manifestHash)
  {
    if (!directory || run<=0 || segment<0 || inputHash.size()!=64 || manifestHash.size()!=64) return false;
    m_directory=directory; m_run=run; m_segment=segment;
    m_inputHash=inputHash; m_manifestHash=manifestHash;
    TDirectory::TContext context(directory);
    m_tree=new TTree("RJTriggerScalersV1","Source-bound cumulative GL1 scaler snapshots; NOT luminosity");
    m_tree->SetDirectory(directory);
    m_tree->Branch("snapshot_id",&m_id,"snapshot_id/l");
    m_tree->Branch("run",&m_run,"run/I");
    m_tree->Branch("segment",&m_segment,"segment/I");
    m_tree->Branch("input_uri_hash",&m_inputHash);
    m_tree->Branch("source_manifest_sha256",&m_manifestHash);
    m_tree->Branch("first_event_sequence",&m_firstEvent,"first_event_sequence/L");
    m_tree->Branch("last_event_sequence",&m_lastEvent,"last_event_sequence/L");
    m_tree->Branch("first_physical_event",&m_firstPhysical,"first_physical_event/L");
    m_tree->Branch("last_physical_event",&m_lastPhysical,"last_physical_event/L");
    m_tree->Branch("first_bco",&m_firstBco,"first_bco/l");
    m_tree->Branch("last_bco",&m_lastBco,"last_bco/l");
    m_tree->Branch("observed_events",&m_observations,"observed_events/l");
    m_tree->Branch("state",&m_state,"state/I");
    m_tree->Branch("packet_version",&m_version,"packet_version/I");
    m_tree->Branch("packet_status",&m_packetStatus,"packet_status/l");
    m_tree->Branch("discontinuity_flags",&m_flags,"discontinuity_flags/i");
    m_tree->Branch("raw",m_raw.data(),"raw[64]/l");
    m_tree->Branch("live",m_live.data(),"live[64]/l");
    m_tree->Branch("scaled",m_scaled.data(),"scaled[64]/l");
    return true;
  }

  bool observe(const RJGl1ScalerContract::Snapshot& s,std::int64_t eventSequence,
               std::uint64_t& id)
  {
    if (!m_tree || m_finished || (m_have && eventSequence<=m_lastEvent)) return false;
    const unsigned flags=RJGl1ScalerContract::discontinuity(m_have ? &m_previous : nullptr,s);
    if (!m_have || !RJGl1ScalerContract::identicalCounters(m_previous,s) || flags!=0)
    {
      if (m_have && !flush()) return false;
      ++m_id; m_firstEvent=eventSequence; m_firstPhysical=s.physical_event;
      m_firstBco=s.bco; m_observations=0; m_flags=0;
      for (int i=0;i<64;++i) { m_raw[i]=s.raw[i]; m_live[i]=s.live[i]; m_scaled[i]=s.scaled[i]; }
      m_state=s.state; m_version=s.version; m_packetStatus=s.packet_status;
    }
    m_lastEvent=eventSequence; m_lastPhysical=s.physical_event; m_lastBco=s.bco;
    ++m_observations; m_flags|=flags; m_previous=s; m_have=true; id=m_id;
    return true;
  }

  bool finish()
  {
    if (!m_tree || m_finished) return false;
    TDirectory::TContext context(m_directory);
    if (m_have && !flush()) return false;
    if (m_tree->Write("RJTriggerScalersV1",TObject::kOverwrite)<=0) return false;
    TNamed contract("rj_trigger_scalers_contract",
      "v1; uint64 raw/live/scaled[64]; lossless counter-state RLE before retention; "
      "join RJEventV1.trigger_scaler_snapshot_id within exact source; "
      "state 0 valid,1 missing,2 unsupported,3 read error,4 packet status bad; "
      "flags 1 counter decrease,2 BCO regression,4 nonincreasing physical event,"
      "8 missing BCO,16 missing physical event; "
      "CUMULATIVE_COUNTERS_NOT_LUMINOSITY; no automatic reset/wrap correction");
    if (contract.Write(contract.GetName(),TObject::kOverwrite)<=0) return false;
    m_finished=true; return true;
  }

 private:
  bool flush() { return m_tree->Fill()>=0; }
  TDirectory* m_directory=nullptr;
  TTree* m_tree=nullptr;
  RJGl1ScalerContract::Snapshot m_previous;
  bool m_have=false,m_finished=false;
  int m_run=0,m_segment=0,m_state=0,m_version=0;
  std::string m_inputHash,m_manifestHash;
  ULong64_t m_id=0,m_observations=0,m_firstBco=0,m_lastBco=0,m_packetStatus=0;
  Long64_t m_firstEvent=0,m_lastEvent=0,m_firstPhysical=0,m_lastPhysical=0;
  UInt_t m_flags=0;
  std::array<ULong64_t,64> m_raw{},m_live{},m_scaled{};
};
}
#endif
