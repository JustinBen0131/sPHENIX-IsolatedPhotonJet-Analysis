#ifndef RJ_CEMC_CALIBRATION_AUDIT_H
#define RJ_CEMC_CALIBRATION_AUDIT_H
#include "RJCemcTowerStatusGuard.h"
#include <caloreco/CaloTowerCalib.h>
#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>
#include <TMD5.h>
#include <cmath>
#include <limits>
#include <nlohmann/json.hpp>

namespace rj_cemc_status {
// Used only for the explicitly bound AuAu DATA reconstruction path. One
// initialization record per job; no per-tower records or extra event output.
inline std::string payload_md5(const std::string& path) {
  std::unique_ptr<TMD5> hash(TMD5::FileChecksum(path.c_str()));
  require(bool(hash), "payload checksum failed: " + path);
  return hash->AsString();
}
inline unsigned validate_zs(const std::string& path, TowerInfoContainer* towers) {
  require(!path.empty() && towers && towers->size()==channel_count, "missing ZS payload or towers");
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  require(file && !file->IsZombie(), "unreadable ZS payload");
  file->Close();
  CDBTTree payload(path);
  unsigned zero_sentinels=0;
  for(unsigned i=0;i<towers->size();++i) {
    const float ratio=payload.GetFloatValue(towers->encode_key(i), "ratio", 0);
    // CaloTowerCalib explicitly interprets a stored zero as unity, not as a
    // missing field. Missing fields are NaN. Preserve that native convention.
    require(std::isfinite(ratio) && ratio>=0, "missing/negative ZS ratio at channel " + std::to_string(i)+" value="+std::to_string(ratio));
    zero_sentinels += ratio==0;
  }
  return zero_sentinels;
}
class CalibrationAudit final : public CaloTowerCalib {
 public:
  explicit CalibrationAudit(std::shared_ptr<State> state)
    : CaloTowerCalib("CEMCCALIB"), state_(std::move(state)) {
    set_detector_type(CaloTowerDefs::CEMC);
    set_doZScrosscalib(true);
    set_doAbortNoZSCalib(true);
  }
  int InitRun(PHCompositeNode* top) override {
    try {
      require(state_ && state_->initialized, "ZS audit requires initialized status guard");
      auto rc=recoConsts::instance();
      require(rc->get_IntFlag("RUNNUMBER")==state_->expected_run &&
              rc->get_uint64Flag("TIMESTAMP")==static_cast<uint64_t>(state_->expected_run), "calibration run/timestamp mismatch");
      const std::string url=CDBInterface::instance()->getUrl("CEMC_ZSCrossCalib");
      auto towers=findNode::getClass<TowerInfoContainer>(top,"TOWERS_CEMC");
      const unsigned zero_sentinels=validate_zs(url,towers);
      // Pin exactly the bytes whose content was checked; no second lookup or
      // silent disable is permitted inside the underlying calibrator.
      set_directURL_ZScrosscalib(url);
      const int code=CaloTowerCalib::InitRun(top);
      require(code==Fun4AllReturnCodes::EVENT_OK,"CaloTowerCalib InitRun failed");
      initialized_=true;
      nlohmann::json receipt={{"schema","CEMCCalibrationInitializationV1"},{"status","PASS"},
        {"run",state_->expected_run},{"globaltag",rc->get_StringFlag("CDB_GLOBALTAG")},
        {"original_input",state_->original_calo_input},{"bad_tower_payload",state_->selected_payload},
        {"bad_tower_md5",payload_md5(state_->selected_payload)},{"zs_payload",url},
        {"zs_md5",payload_md5(url)},{"zs_channels_validated",channel_count},{"zs_enabled",true},
        {"zs_zero_means_unity_channels",zero_sentinels},
        {"status_action",state_->restoring_missing_status?"RESTORE_MISSING_HOT_ONLY":"PRESERVE_VALIDATE"}};
      std::cout<<"CEMC_CALIBRATION_INIT\t"<<receipt.dump()<<std::endl;
      return code;
    } catch(const std::exception& e) {
      std::cerr<<"CEMC_CALIBRATION_FATAL "<<e.what()<<std::endl;
      gSystem->Exit(1);return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  int process_event(PHCompositeNode* top) override {
    if(!initialized_) { gSystem->Exit(1);return Fun4AllReturnCodes::ABORTRUN; }
    const int code=CaloTowerCalib::process_event(top);
    if(code!=Fun4AllReturnCodes::EVENT_OK) { gSystem->Exit(1);return code; }
    ++events_;return code;
  }
  int End(PHCompositeNode* top) override {
    const int code=CaloTowerCalib::End(top);
    std::cout<<"CEMC_CALIBRATION_END\t"<<nlohmann::json({{"run",state_->expected_run},
      {"initialized",initialized_},{"events",events_},{"returncode",code}}).dump()<<std::endl;
    return code;
  }
 private:
  std::shared_ptr<State> state_;
  bool initialized_=false;
  unsigned long long events_=0;
};
}
#endif
