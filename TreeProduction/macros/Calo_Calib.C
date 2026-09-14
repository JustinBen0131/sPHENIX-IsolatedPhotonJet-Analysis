#ifndef CALO_CALIB_H
#define CALO_CALIB_H

#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterDeadHotMask.h>
#include <caloreco/RawClusterPositionCorrection.h>
#include "calo/RJCemcTowerStatusGuard.h"
#include "calo/RJCemcCalibrationAudit.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <calostatusskimmer/CaloStatusSkimmer.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>  // for Fun4AllServer
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/RunnumberRange.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <TSystem.h>  // for gSystem

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <string>
#include <vector>

R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4allutils.so)
R__LOAD_LIBRARY(libCaloStatusSkimmer.so)

namespace
{
  bool rj_env_truthy(const char* name)
  {
    const char* raw = getenv(name);
    if (!raw)
    {
      return false;
    }
    std::string flag(raw);
    std::transform(flag.begin(), flag.end(), flag.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return (flag == "1" || flag == "true" || flag == "yes" || flag == "on");
  }

  class RJCaloTowerStatusAudit : public SubsysReco
  {
   public:
    explicit RJCaloTowerStatusAudit(const std::string& name,
                                    const std::vector<std::string>& nodes)
      : SubsysReco(name)
      , m_nodes(nodes)
    {
    }

    int process_event(PHCompositeNode* topNode) override
    {
      if (m_done)
      {
        return Fun4AllReturnCodes::EVENT_OK;
      }
      m_done = true;

      std::cout << "[RJCaloTowerStatusAudit] first-event calibrated TowerInfo status audit" << std::endl;
      for (const auto& node : m_nodes)
      {
        TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, node);
        if (!towers)
        {
          std::cout << "[RJCaloTowerStatusAudit] node=" << node << " missing" << std::endl;
          continue;
        }

        unsigned int present = 0;
        unsigned int good = 0;
        unsigned int bad = 0;
        unsigned int hot = 0;
        unsigned int badChi2 = 0;
        unsigned int noCalib = 0;
        unsigned int notInstr = 0;

        const unsigned int nchannels = towers->size();
        for (unsigned int channel = 0; channel < nchannels; ++channel)
        {
          TowerInfo* tower = towers->get_tower_at_channel(channel);
          if (!tower)
          {
            continue;
          }
          ++present;
          if (tower->get_isGood())
          {
            ++good;
          }
          else
          {
            ++bad;
          }
          if (tower->get_isHot())
          {
            ++hot;
          }
          if (tower->get_isBadChi2())
          {
            ++badChi2;
          }
          if (tower->get_isNoCalib())
          {
            ++noCalib;
          }
          if (tower->get_isNotInstr())
          {
            ++notInstr;
          }
        }

        std::cout << "[RJCaloTowerStatusAudit] node=" << node
                  << " channels=" << nchannels
                  << " present=" << present
                  << " good=" << good
                  << " bad=" << bad
                  << " hot=" << hot
                  << " badChi2=" << badChi2
                  << " noCalib=" << noCalib
                  << " notInstr=" << notInstr
                  << std::endl;
      }
      return Fun4AllReturnCodes::EVENT_OK;
    }

   private:
    std::vector<std::string> m_nodes;
    bool m_done = false;
  };
}

void Process_Calo_Calib(std::shared_ptr<rj_cemc_status::State> cemcStatusState = {})
{
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc = recoConsts::instance();

  /////////////////
  // set MC or data
  bool isSim = true;
  int data_sim_runnumber_thres = 1000;
  int runnumber = rc->get_uint64Flag("TIMESTAMP");
  if (rc->get_uint64Flag("TIMESTAMP") > data_sim_runnumber_thres)
  {
    isSim = false;
  }
  std::cout << "Calo Calib uses runnumber " << rc->get_uint64Flag("TIMESTAMP") << std::endl;

  // ------------------------------------------------------------------
  // UPDATE (embedded-only):
  // Detect the embedded-simulation pipeline from RJ_DATASET so we can
  // keep Process_Calo_Calib() in the existing data-like branch
  // (TIMESTAMP > 1000) while applying a very small exception only for
  // the embedded path below.
  //
  // This does NOT affect any non-embedded pipeline.
  // ------------------------------------------------------------------
  bool isSimEmbedded = false;
  if (const char* ds = getenv("RJ_DATASET"))
  {
    std::string dataset(ds);
    if (dataset == "isSimEmbedded" || dataset == "issimembedded" || dataset == "simembedded" ||
        dataset == "isSimEmbeddedInclusive" || dataset == "issimembeddedinclusive" ||
        dataset == "simembeddedinclusive")
    {
      isSimEmbedded = true;
    }
  }

  bool isScaledTriggerStudyOnly = rj_env_truthy("RJ_SCALED_TRIGGER_STUDY_ONLY");
  bool skipCaloStatusSkimmer = rj_env_truthy("RJ_SKIP_CALO_STATUS_SKIMMER");
  bool skipLegacyCaloTowerStatus = isScaledTriggerStudyOnly || rj_env_truthy("RJ_SKIP_CALO_TOWER_STATUS");
  bool forceEmbeddedCaloTowerStatus = rj_env_truthy("RJ_FORCE_CALO_TOWER_STATUS_FOR_EMBEDDED");
  bool auditCalibTowerStatus = rj_env_truthy("RJ_CALO_STATUS_AUDIT");
  const char* statusInputPrefixRaw = getenv("RJ_CALO_TOWER_STATUS_INPUT_PREFIX");
  const std::string statusInputPrefix = statusInputPrefixRaw ? std::string(statusInputPrefixRaw) : std::string();
  const bool useStatusInputPrefix = !statusInputPrefix.empty();

  ///////////////////////////////////////////////
  // Remove incomplete events from event combiner
  if (!isSim && !isSimEmbedded && !skipCaloStatusSkimmer)
  {
    CaloStatusSkimmer *css = new CaloStatusSkimmer("CaloStatusSkimmer");
    se->registerSubsystem(css);
  }
  else if (isSimEmbedded)
  {
    std::cout << "[Process_Calo_Calib][isSimEmbedded] skipping CaloStatusSkimmer "
                 "(data event-combiner completeness guard is not used for embedded SIM)" << std::endl;
  }
  else if (skipCaloStatusSkimmer)
  {
    std::cout << "[Process_Calo_Calib] skipping CaloStatusSkimmer "
                 "(CALOFITTING/TowerInfo input is analyzed as the primary data stream)" << std::endl;
  }

  //////////////////////
  // Input geometry node
  std::cout << "Adding Geometry file" << std::endl;
  Fun4AllInputManager *ingeo = new Fun4AllRunNodeInputManager("DST_GEO");
  std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
  ingeo->AddFile(geoLocation);
  se->registerInputManager(ingeo);

  CaloTowerDefs::BuilderType buildertype = CaloTowerDefs::kPRDFTowerv4;

  // build ZDC towers
  if (isSimEmbedded)
  {
    std::cout << "[Process_Calo_Calib][isSimEmbedded] skipping ZDCBUILDER "
                 "(ZDC/minimum-bias calibration is data-only for this embedded SIM path)" << std::endl;
  }
  else
  {
    CaloTowerBuilder *caZDC = new CaloTowerBuilder("ZDCBUILDER");
    caZDC->set_detector_type(CaloTowerDefs::ZDC);
    caZDC->set_builder_type(buildertype);
    if ((runnumber > RunnumberRange::RUN2PP_FIRST && runnumber < RunnumberRange::RUN2PP_LAST) || (runnumber > RunnumberRange::RUN3PP_FIRST && runnumber < RunnumberRange::RUN3PP_LAST))
    {
      caZDC->set_processing_type(CaloWaveformProcessing::FAST);
    }
    else
    {
      caZDC->set_processing_type(CaloWaveformProcessing::FUNCFIT);
      caZDC->set_funcfit_type(2);
    }
    caZDC->set_nsamples(16);
    caZDC->set_offlineflag();
    se->registerSubsystem(caZDC);
  }

  //////////////////////////////
  // set statuses on raw towers
  // AuAu paired-DATA has raw TowerInfo status copied by CaloTowerCalib.
  // Validate (or explicitly recover) that status before any calo consumer.
  // Do not also run a resetting legacy status setter under a different tag.
  if (cemcStatusState)
  {
    if (isSim || isSimEmbedded || !skipLegacyCaloTowerStatus || useStatusInputPrefix)
    {
      std::cerr << "CEMC_STATUS_CONTRACT_FATAL unsupported paired-DATA status configuration" << std::endl;
      gSystem->Exit(1);
    }
    se->registerSubsystem(new rj_cemc_status::Guard(cemcStatusState, false));
  }

  if (isSimEmbedded && !forceEmbeddedCaloTowerStatus)
  {
    std::cout << "[Process_Calo_Calib][isSimEmbedded] skipping CaloTowerStatus setters "
                 "(embedded SIM path uses producer tower-quality state)" << std::endl;
  }
  else if (skipLegacyCaloTowerStatus)
  {
    std::cout << "[Process_Calo_Calib] skipping CaloTowerStatus setters "
                 "(CALOFITTING/TowerInfo input has no legacy TOWERS_ node)" << std::endl;
  }
  else
  {
    std::cout << "status setters" << std::endl;
    CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
    statusEMC->set_detector_type(CaloTowerDefs::CEMC);
    if (useStatusInputPrefix)
    {
      statusEMC->set_inputNodePrefix(statusInputPrefix);
      std::cout << "[Process_Calo_Calib] CEMCSTATUS input prefix: "
                << statusInputPrefix << std::endl;
    }
    // MC Towers Status
    if (isSim)
    {
      // Uses threshold of 50% for towers be considered frequently bad.
      std::string calibName_hotMap = "CEMC_hotTowers_status";
      /* Systematic options (to be used as needed). */
      /* Uses threshold of 40% for towers be considered frequently bad. */
      // std::string calibName_hotMap = "CEMC_hotTowers_status_40";

      /* Uses threshold of 60% for towers be considered frequently bad. */
      // std::string calibName_hotMap = "CEMC_hotTowers_status_60";

      std::string calibdir = CDBInterface::instance()->getUrl(calibName_hotMap);
      statusEMC->set_directURL_hotMap(calibdir);
    }
    if (auditCalibTowerStatus)
    {
      statusEMC->Verbosity(1);
    }
    se->registerSubsystem(statusEMC);

    CaloTowerStatus *statusHCalIn = new CaloTowerStatus("HCALINSTATUS");
    statusHCalIn->set_detector_type(CaloTowerDefs::HCALIN);
    if (useStatusInputPrefix)
    {
      statusHCalIn->set_inputNodePrefix(statusInputPrefix);
      std::cout << "[Process_Calo_Calib] HCALINSTATUS input prefix: "
                << statusInputPrefix << std::endl;
    }
    if (auditCalibTowerStatus)
    {
      statusHCalIn->Verbosity(1);
    }
    se->registerSubsystem(statusHCalIn);

    CaloTowerStatus *statusHCALOUT = new CaloTowerStatus("HCALOUTSTATUS");
    statusHCALOUT->set_detector_type(CaloTowerDefs::HCALOUT);
    if (useStatusInputPrefix)
    {
      statusHCALOUT->set_inputNodePrefix(statusInputPrefix);
      std::cout << "[Process_Calo_Calib] HCALOUTSTATUS input prefix: "
                << statusInputPrefix << std::endl;
    }
    if (auditCalibTowerStatus)
    {
      statusHCALOUT->Verbosity(1);
    }
    se->registerSubsystem(statusHCALOUT);
  }

  ////////////////////
  // Calibrate towers
  std::cout << "Calibrating EMCal" << std::endl;
  CaloTowerCalib *calibEMC = cemcStatusState
    ? static_cast<CaloTowerCalib*>(new rj_cemc_status::CalibrationAudit(cemcStatusState))
    : new CaloTowerCalib("CEMCCALIB");
  calibEMC->set_detector_type(CaloTowerDefs::CEMC);
  if (useStatusInputPrefix)
  {
    calibEMC->set_inputNodePrefix(statusInputPrefix);
  }
  se->registerSubsystem(calibEMC);

  std::cout << "Calibrating OHcal" << std::endl;
  CaloTowerCalib *calibOHCal = new CaloTowerCalib("HCALOUT");
  calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
  if (useStatusInputPrefix)
  {
    calibOHCal->set_inputNodePrefix(statusInputPrefix);
  }
  se->registerSubsystem(calibOHCal);

  std::cout << "Calibrating IHcal" << std::endl;
  CaloTowerCalib *calibIHCal = new CaloTowerCalib("HCALIN");
  calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
  if (useStatusInputPrefix)
  {
    calibIHCal->set_inputNodePrefix(statusInputPrefix);
  }
  se->registerSubsystem(calibIHCal);

  if (cemcStatusState)
  {
    se->registerSubsystem(new rj_cemc_status::Guard(cemcStatusState, true));
  }

  if ((!isSimEmbedded || forceEmbeddedCaloTowerStatus) && auditCalibTowerStatus)
  {
    std::cout << "[Process_Calo_Calib] auditing TOWERINFO_CALIB_* status after CaloTowerCalib copy "
                 "and before RawClusterBuilderTemplate" << std::endl;
    se->registerSubsystem(new RJCaloTowerStatusAudit(
        "RJCaloTowerStatusAudit",
        {"TOWERINFO_CALIB_CEMC", "TOWERINFO_CALIB_HCALIN", "TOWERINFO_CALIB_HCALOUT"}));
  }

  ////////////////
  // MC Calibration
  if (isSim && rc->get_uint64Flag("TIMESTAMP") < 28)  // in run28 and beyond we moved the MC calibration into the waveformsim module for data embedding
  {
    std::string MC_Calib = CDBInterface::instance()->getUrl("CEMC_MC_RECALIB");
    if (MC_Calib.empty())
    {
      std::cout << "No MC calibration found :( )" << std::endl;
      gSystem->Exit(0);
    }
    CaloTowerCalib *calibEMC_MC = new CaloTowerCalib("CEMCCALIB_MC");
    calibEMC_MC->set_detector_type(CaloTowerDefs::CEMC);
    calibEMC_MC->set_inputNodePrefix("TOWERINFO_CALIB_");
    calibEMC_MC->set_outputNodePrefix("TOWERINFO_CALIB_");
    calibEMC_MC->set_directURL(MC_Calib);
    calibEMC_MC->set_doCalibOnly(true);
    se->registerSubsystem(calibEMC_MC);
  }

  //////////////////
  // Clusters
  std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->set_threshold_energy(0.070);  // for when using basic calibration
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->set_UseTowerInfo(1);   // to use towerinfo objects rather than old RawTower
  ClusterBuilder->set_UseAltZVertex(1);  // Use MBD Vertex for vertex-based corrections
  se->registerSubsystem(ClusterBuilder);
}

#endif
