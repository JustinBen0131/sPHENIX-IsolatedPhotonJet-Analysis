#ifndef RJ_TRIGGER_RUN_INFO_V1_H
#define RJ_TRIGGER_RUN_INFO_V1_H
#include <TDirectory.h>
#include <TNamed.h>
#include <TTree.h>
#include <cstdint>
#include <limits>
#include <string>

namespace RJTriggerRunInfoV1
{
// One run-configuration table per physical source, not 64 names/prescales
// repeated on every event. Run-total counters are distinct from event GL1
// snapshots. Initial hardware prescale and run-averaged prescale stay separate.
class Writer
{
 public:
  bool initialize(TDirectory* dir,int run,int segment,const std::string& input,
                  const std::string& manifest)
  {
    if(!dir || run<=0 || segment<0 || input.size()!=64 || manifest.size()!=64) return false;
    TDirectory::TContext context(dir);
    m_dir=dir; m_run=run; m_segment=segment; m_input=input; m_manifest=manifest;
    m_tree=new TTree("RJTriggerRunInfoV1","Native run trigger configuration; source scoped");
    m_tree->SetDirectory(dir);
    m_tree->Branch("run",&m_run,"run/I");
    m_tree->Branch("segment",&m_segment,"segment/I");
    m_tree->Branch("input_uri_hash",&m_input);
    m_tree->Branch("source_manifest_sha256",&m_manifest);
    m_tree->Branch("bit",&m_bit,"bit/I");
    m_tree->Branch("name",&m_name);
    m_tree->Branch("state",&m_state,"state/I");
    m_tree->Branch("initial_prescale",&m_initial,"initial_prescale/I");
    m_tree->Branch("prescale",&m_prescale,"prescale/D");
    m_tree->Branch("run_raw",&m_raw,"run_raw/l");
    m_tree->Branch("run_live",&m_live,"run_live/l");
    m_tree->Branch("run_scaled",&m_scaled,"run_scaled/l");
    return true;
  }

  template<class RunInfo> bool observe(const RunInfo* info)
  {
    if(!m_tree) return false;
    if(m_written) return true; // RUN-node configuration is constant per source
    const bool supported=info && std::string(info->ClassName())=="TriggerRunInfov1";
    TDirectory::TContext context(m_dir);
    for(m_bit=0;m_bit<64;++m_bit)
    {
      m_name=""; m_initial=-1; m_prescale=std::numeric_limits<double>::quiet_NaN();
      m_raw=m_live=m_scaled=0;
      m_state=info ? 2 : 1; // 0 known named bit,1 missing node,2 unsupported,3 unassigned
      if(supported)
      {
        m_name=info->getTriggerName(m_bit);
        m_initial=info->getInitialPrescaleByBit(m_bit);
        m_prescale=info->getPrescaleByBit(m_bit);
        m_raw=info->getRawScalersByBit(m_bit);
        m_live=info->getLiveScalersByBit(m_bit);
        m_scaled=info->getScalersByBit(m_bit);
        m_state=(m_name.empty() || m_name=="unknown") ? 3 : 0;
      }
      if(m_tree->Fill()<0) return false;
    }
    m_written=true;
    return true;
  }

  bool finish()
  {
    if(!m_tree) return false;
    TDirectory::TContext context(m_dir);
    if(m_tree->Write("RJTriggerRunInfoV1",TObject::kOverwrite)<=0) return false;
    TNamed contract("rj_trigger_retention_contract",
      "direct_gl1_v1; live/scaled event decisions and independent raw/live/scaled counters; "
      "TriggerRunInfo names, initial hardware prescale, native run-averaged prescale and run totals; "
      "no DATA-trigger cut on simulation; workfest capture retains every producer encounter; "
      "object cuts and upstream DST/reconstruction population remain separately defined; "
      "TriggerVector is a legacy live alias in Gl1Packetv2/v3, not raw bits; "
      "counter differences require source/run coverage and discontinuity reconciliation");
    return contract.Write(contract.GetName(),TObject::kOverwrite)>0;
  }
 private:
  TDirectory* m_dir=nullptr;
  TTree* m_tree=nullptr;
  bool m_written=false;
  int m_run=0,m_segment=0,m_bit=0,m_state=1,m_initial=-1;
  std::string m_input,m_manifest,m_name;
  double m_prescale=0;
  ULong64_t m_raw=0,m_live=0,m_scaled=0;
};
}
#endif
