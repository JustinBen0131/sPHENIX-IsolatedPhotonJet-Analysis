//======================================================================
//  Fun4All_recoilJets_unified_impl.C
//  --------------------------------------------------------------------
#pragma once
#if !defined(RJ_UNIFIED_ANALYSIS_PP) && !defined(RJ_UNIFIED_ANALYSIS_AUAU)
  #error "Define RJ_UNIFIED_ANALYSIS_PP or RJ_UNIFIED_ANALYSIS_AUAU before including Fun4All_recoilJets_unified_impl.C"
#endif
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

//–––– Standard Fun4All ––––––––––––––––––––––––––––––––––––––
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <jetbase/JetContainer.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <mbd/MbdPmtContainer.h>
#include <ffarawobjects/Gl1Packet.h>
#include <ffaobjects/EventHeader.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/CaloTowerBuilder.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <frog/FROG.h>
#include <calotrigger/MinimumBiasClassifier.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/CDBInterface.h>
#include <fun4allutils/TimerStats.h>
#include <clusteriso/ClusterIso.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <caloreco/CaloGeomMapping.h>
#include <caloreco/RawClusterPositionCorrection.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterBuilderTopo.h>
#include <calobase/RawTowerGeom.h>
#include <caloreco/RawTowerCalibration.h>
#include <calowaveformsim/CaloWaveformSim.h>
#include <caloreco/PhotonClusterBuilder.h>
#include <jetbase/Jet.h>
#include <g4jets/TruthJetInput.h>

#include <jetbase/FastJetOptions.h>
#include <caloana/RJFastJetAlgoSubArea.h>
#include <caloana/RJReplayFoundationV1.h>
#include <globalvertex/GlobalVertexReco.h>
#include <caloreco/CaloTowerCalib.h>

#include <centrality/CentralityReco.h>
#include <centrality/CentralityInfo.h>
#include <g4mbd/MbdDigitization.h>
#include <mbd/MbdEvent.h>
#include <mbd/MbdReco.h>
#include <zdcinfo/ZdcReco.h>
#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>
#include <jetbase/FastJetOptions.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbase/JetCalib.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/TowerBackground.h>
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
#include <eventplaneinfo/EventPlaneReco.h>
#include <caloana/RecoilJets_AuAu.h>
#else
#include <caloana/RecoilJets.h>
#endif

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>     // getenv
#include <algorithm>   // std::transform
#include <cctype>      // std::tolower
#include <iomanip>     // std::setw, std::setprecision
#include <dlfcn.h>     // dlopen RTLD_NOLOAD
#include <typeinfo>    // typeid RTTI probe
#include <TSystem.h>   // gSystem, GetBuildArch/Compiler info
#include <csignal>     // signal handlers (debug backtrace)
#include <cerrno>
#include <execinfo.h>  // backtrace
#include <unistd.h>    // STDERR_FILENO
#include <cstdio>      // snprintf
#include <limits>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TNamed.h>
#include <TObject.h>
#include <TRandom3.h>
#include "Calo_Calib.C"

// Calo_Calib.C loads the private CaloReco stack once.  Do not explicitly
// reload libcalo_reco/libcalo_io here: ROOT will re-register the CaloBase
// dictionaries and pp-SIM G4 rebuild jobs can abort before event processing.

#if defined(__has_include)
#  if __has_include(<GlobalVariables.C>) && __has_include(<G4_Input.C>) && __has_include(<G4_CEmc_Spacal.C>) && __has_include(<G4_HcalIn_ref.C>) && __has_include(<G4_HcalOut_ref.C>) && __has_include(<G4_Mbd.C>) && __has_include(<G4_RunSettings.C>)
#    include <GlobalVariables.C>
#    include <G4_Input.C>
#    include <G4_CEmc_Spacal.C>
#    include <G4_HcalIn_ref.C>
#    include <G4_HcalOut_ref.C>
#    include <G4_Mbd.C>
#    include <G4_RunSettings.C>
#    define RJ_HAS_SPHENIX_G4_INPUT_MACROS 1
#  endif
#endif
#ifndef RJ_HAS_SPHENIX_G4_INPUT_MACROS
#  define RJ_HAS_SPHENIX_G4_INPUT_MACROS 0
#endif

// CaloReco/CaloIO are intentionally provided by Calo_Calib.C above.
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
R__LOAD_LIBRARY(libRecoilJetsAuAu.so)
#else
R__LOAD_LIBRARY(libRecoilJets.so)
#endif
R__LOAD_LIBRARY(libclusteriso.so)
R__LOAD_LIBRARY(libjetbase.so)


// Then load the rest of the environment stack
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libglobalvertex.so)
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
R__LOAD_LIBRARY(libeventplaneinfo.so)
#endif
R__LOAD_LIBRARY(libcentrality.so)      // always
R__LOAD_LIBRARY(libcentrality_io.so)   // if you instantiate CentralityReco
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libzdcinfo.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libg4mbd.so)

//======================================================================
//  Convenience helpers
//======================================================================
namespace detail
{
  // BEGIN RJ_FUN4ALL_TERMINAL_CLASSIFIER_V1_PURE
  struct Fun4AllTerminalStatusWitnessV1
  {
    int nEventsRequested = 0;
    int runRc = 0;
    int endRc = 0;
    int eventOk = 0;
    int abortProcessingCount = 0;
    int abortRunCount = 0;
    int registeredInputManagers = 0;
    int ordinaryInputManagers = 0;
    int exhaustedOrdinaryInputManagers = 0;
    int openOrdinaryInputManagers = 0;
    int nonemptyOrdinaryFileLists = 0;
    int permittedRepeatingManagerExpected = 0;
    int permittedRepeatingManagerMatches = 0;
    int permittedRunNodeInputManagers = 0;
  };

  enum class Fun4AllTerminalStatusClassV1
  {
    kEventOk,
    kVerifiedMultiInputEof,
    kFail
  };

  struct Fun4AllTerminalStatusDecisionV1
  {
    Fun4AllTerminalStatusClassV1 classification =
      Fun4AllTerminalStatusClassV1::kFail;
    const char* status = "FAIL";
    const char* reason = "unclassified";
  };

  inline Fun4AllTerminalStatusDecisionV1 classify_fun4all_terminal_status(
      const Fun4AllTerminalStatusWitnessV1& witness)
  {
    if (witness.runRc == witness.eventOk &&
        witness.endRc == witness.eventOk)
    {
      return {
        Fun4AllTerminalStatusClassV1::kEventOk,
        "PASS",
        "EVENT_OK"};
    }
    if (witness.endRc != witness.eventOk)
    {
      return {
        Fun4AllTerminalStatusClassV1::kFail,
        "FAIL",
        "END_NONZERO"};
    }
    if (witness.abortProcessingCount != 0 ||
        witness.abortRunCount != 0)
    {
      return {
        Fun4AllTerminalStatusClassV1::kFail,
        "FAIL",
        "ABORT_STATISTIC_NONZERO"};
    }
    if (witness.permittedRepeatingManagerExpected < 0 ||
        witness.permittedRepeatingManagerExpected > 1 ||
        witness.permittedRepeatingManagerMatches !=
          witness.permittedRepeatingManagerExpected)
    {
      return {
        Fun4AllTerminalStatusClassV1::kFail,
        "FAIL",
        "REPEATING_MANAGER_POINTER_MISMATCH"};
    }
    if (witness.ordinaryInputManagers <= 0 ||
        witness.permittedRunNodeInputManagers < 0 ||
        witness.registeredInputManagers !=
          witness.ordinaryInputManagers +
            witness.permittedRepeatingManagerExpected +
            witness.permittedRunNodeInputManagers)
    {
      return {
        Fun4AllTerminalStatusClassV1::kFail,
        "FAIL",
        "INPUT_MANAGER_ACCOUNTING_MISMATCH"};
    }
    if (witness.exhaustedOrdinaryInputManagers !=
          witness.ordinaryInputManagers ||
        witness.openOrdinaryInputManagers != 0 ||
        witness.nonemptyOrdinaryFileLists != 0)
    {
      return {
        Fun4AllTerminalStatusClassV1::kFail,
        "FAIL",
        "ORDINARY_INPUT_NOT_EXHAUSTED"};
    }
    if (witness.runRc != -witness.ordinaryInputManagers)
    {
      return {
        Fun4AllTerminalStatusClassV1::kFail,
        "FAIL",
        "EOF_SUM_MISMATCH"};
    }
    return {
      Fun4AllTerminalStatusClassV1::kVerifiedMultiInputEof,
      "PASS",
      witness.nEventsRequested == 0
        ? "VERIFIED_MULTI_INPUT_EOF"
        : "VERIFIED_MULTI_INPUT_EOF_BEFORE_UPPER_BOUND"};
  }
  // END RJ_FUN4ALL_TERMINAL_CLASSIFIER_V1_PURE

  /// Throw a nicely formatted exception on unrecoverable error
  [[noreturn]] void bail(const std::string& msg)
  {
    std::ostringstream oss;
    oss << "\n[FATAL] Fun4All_recoilJets :: " << msg << '\n';
    throw std::runtime_error(oss.str());
  }

  /// Preserve the two terminal Fun4All return codes even in quiet Condor
  /// mode. ScopedSilence redirects the C++ iostream buffers only; the narrow
  /// C stderr record remains visible without restoring ordinary batch chatter.
  inline void enforce_fun4all_status(const char* path,
                                     Fun4AllServer* server,
                                     const int nEvents,
                                     const int runRc,
                                     const int endRc,
                                     const Fun4AllInputManager*
                                       permittedRepeatingManager)
  {
    Fun4AllTerminalStatusWitnessV1 witness;
    witness.nEventsRequested = nEvents;
    witness.runRc = runRc;
    witness.endRc = endRc;
    witness.eventOk = Fun4AllReturnCodes::EVENT_OK;
    witness.permittedRepeatingManagerExpected =
      permittedRepeatingManager ? 1 : 0;

    if (server)
    {
      witness.abortProcessingCount =
        server->retcodestats(Fun4AllReturnCodes::ABORTPROCESSING);
      witness.abortRunCount =
        server->retcodestats(Fun4AllReturnCodes::ABORTRUN);
      auto* syncManager = server->getSyncManager();
      if (syncManager)
      {
        const auto& inputManagers = syncManager->GetInputManagers();
        witness.registeredInputManagers =
          static_cast<int>(inputManagers.size());
        int inputManagerIndex = 0;
        for (const auto* inputManager : inputManagers)
        {
          const bool isPermittedRepeating =
            inputManager == permittedRepeatingManager;
          const bool isPermittedRunNode =
            dynamic_cast<const Fun4AllRunNodeInputManager*>(inputManager) !=
              nullptr;
          const bool isOpen = inputManager && inputManager->IsOpen();
          const bool fileListEmpty =
            inputManager && inputManager->FileListEmpty();
          std::fprintf(
            stderr,
            "RECOILJETS_FUN4ALL_INPUT_STATUS_V1"
            " path=%s index=%d name=%s repeating=%d run_node_auxiliary=%d"
            " open=%d file_list_empty=%d\n",
            path,
            inputManagerIndex,
            inputManager ? inputManager->Name().c_str() : "<null>",
            isPermittedRepeating ? 1 : 0,
            isPermittedRunNode ? 1 : 0,
            isOpen ? 1 : 0,
            fileListEmpty ? 1 : 0);
          ++inputManagerIndex;
          if (inputManager == permittedRepeatingManager)
          {
            ++witness.permittedRepeatingManagerMatches;
            continue;
          }
          if (isPermittedRunNode)
          {
            ++witness.permittedRunNodeInputManagers;
            continue;
          }
          ++witness.ordinaryInputManagers;
          if (!inputManager)
          {
            continue;
          }
          if (isOpen)
          {
            ++witness.openOrdinaryInputManagers;
          }
          if (!fileListEmpty)
          {
            ++witness.nonemptyOrdinaryFileLists;
          }
          if (!isOpen && fileListEmpty)
          {
            ++witness.exhaustedOrdinaryInputManagers;
          }
        }
      }
    }

    const auto decision = classify_fun4all_terminal_status(witness);
    const bool ok =
      decision.classification != Fun4AllTerminalStatusClassV1::kFail;
    std::fprintf(
      stderr,
      "RECOILJETS_FUN4ALL_STATUS_V2 path=%s run_rc=%d end_rc=%d"
      " n_events=%d registered_inputs=%d ordinary_inputs=%d"
      " exhausted_ordinary_inputs=%d open_ordinary_inputs=%d"
      " nonempty_ordinary_file_lists=%d"
      " repeating_expected=%d repeating_matches=%d"
      " run_node_auxiliary_inputs=%d"
      " abort_processing=%d abort_run=%d status=%s reason=%s\n",
      path,
      runRc,
      endRc,
      nEvents,
      witness.registeredInputManagers,
      witness.ordinaryInputManagers,
      witness.exhaustedOrdinaryInputManagers,
      witness.openOrdinaryInputManagers,
      witness.nonemptyOrdinaryFileLists,
      witness.permittedRepeatingManagerExpected,
      witness.permittedRepeatingManagerMatches,
      witness.permittedRunNodeInputManagers,
      witness.abortProcessingCount,
      witness.abortRunCount,
      decision.status,
      decision.reason);
    std::fflush(stderr);
    if (ok) return;

    std::fprintf(
      stderr,
      "[FATAL] Fun4All_recoilJets status failure:"
      " path=%s run_rc=%d end_rc=%d reason=%s\n",
      path,
      runRc,
      endRc,
      decision.reason);
    std::fflush(stderr);
    if (gSystem) gSystem->Exit(90);
    throw std::runtime_error("Fun4All returned a nonzero terminal status");
  }

  /// Trim whitespace from both ends (for robust list-file parsing)
  inline std::string trim(std::string s)
  {
    const char* ws = " \t\r\n";
    s.erase(0, s.find_first_not_of(ws));
    s.erase(s.find_last_not_of(ws) + 1);
    return s;
  }

  inline std::string fmt(double x, int precision = 6)
  {
    std::ostringstream os;
    os << std::fixed << std::setprecision(precision) << x;
    std::string s = os.str();
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    return s.empty() ? "0" : s;
  }

}

namespace detail
{
  inline RJFastJetAlgoSubArea* fjAlgo(const float R, const bool calculateArea = false)
  {
    FastJetOptions o{};               // IMPORTANT: value-initialize ALL fields to safe defaults
    o.algo            = Jet::ANTIKT;  // algorithm
    o.jet_R           = R;            // jet radius
    o.use_jet_min_pt  = true;         // enable a ptmin
    o.jet_min_pt      = 0.0f;         // ptmin value
    o.calc_area       = calculateArea; // genuine FastJet active area when retained
    o.verbosity       = 0;            // quiet
    return new RJFastJetAlgoSubArea(o);
  }
}

namespace rj_centrality_binding
{
  struct Binding
  {
    bool apply_local = false;
    int run = 0;
    std::string divisions, scale, vertex_scale;
    std::string divisions_sha256, scale_sha256, vertex_scale_sha256;
  };

  inline std::string env(const char* key)
  { const char* raw = std::getenv(key); return raw ? detail::trim(raw) : ""; }

  inline void require_run_file(const std::string& path,
                               const std::string& basename)
  {
    if (path.empty() || path.front() != '/' ||
        path.substr(path.find_last_of('/') + 1) != basename)
      detail::bail("centrality payload path is not absolute/run-bound: " + path);
    std::ifstream probe(path, std::ios::in | std::ios::binary);
    if (!probe.good())
      detail::bail("centrality payload unreadable: " + path);
  }

  inline void require_sha256(const std::string& digest)
  {
    if (digest.size() != 64 ||
        !std::all_of(digest.begin(), digest.end(),
                     [](unsigned char c) { return std::isxdigit(c); }))
      detail::bail("centrality SHA256 must be 64 hex characters");
  }

  inline void verify_file_sha256(const std::string& path, std::string expected)
  {
    if (!std::all_of(path.begin(), path.end(), [](unsigned char c) {
          return std::isalnum(c) || c == '/' || c == '.' || c == '_' || c == '-';
        }))
      detail::bail("unsafe centrality path for SHA256");
    std::transform(expected.begin(), expected.end(), expected.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    require_sha256(expected);
#ifdef __APPLE__
    const std::string command = "/usr/bin/shasum -a 256 -- " + path;
#else
    const std::string command = "sha256sum -- " + path;
#endif
    FILE* pipe = ::popen(command.c_str(), "r");
    char output[256] = {};
    if (!pipe || !std::fgets(output, sizeof(output), pipe))
    {
      if (pipe) ::pclose(pipe);
      detail::bail("centrality SHA256 tool failed: " + path);
    }
    errno = 0;
    const int status = ::pclose(pipe);
    const int close_errno = errno;
    std::string observed(output);
    observed = observed.substr(0, observed.find_first_of(" \t\r\n"));
    require_sha256(observed);
    // ROOT/Fun4All can reap the short-lived hash child through SIGCHLD before
    // pclose waits for it. ECHILD is not evidence of a hash mismatch: retain
    // the exact 64-hex digest comparison and reject every other tool error.
    const bool externally_reaped = status == -1 && close_errno == ECHILD;
    if ((status != 0 && !externally_reaped) || observed != expected)
      detail::bail("centrality payload SHA256 mismatch: " + path +
                   " expected=" + expected + " observed=" + observed +
                   " tool_status=" + std::to_string(status) +
                   " close_errno=" + std::to_string(close_errno));
    std::cout << "CENTRALITY_HASH path=" << path
              << " expected_sha256=" << expected
              << " verified_sha256=" << observed
              << " tool_status=" << status << " close_errno=" << close_errno
              << " externally_reaped=" << externally_reaped << std::endl;
  }

  inline Binding resolve(bool is_auau_data, bool centrality_will_run, int run)
  {
    Binding binding;
    binding.run = run;
    std::string source = env("RJ_AUAU_CENTRALITY_SOURCE");
    std::transform(source.begin(), source.end(), source.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    binding.divisions = env("RJ_AUAU_CENTRALITY_DIVS");
    binding.scale = env("RJ_AUAU_CENTRALITY_SCALE");
    binding.vertex_scale = env("RJ_AUAU_CENTRALITY_VERTEX_SCALE");
    binding.divisions_sha256 = env("RJ_AUAU_CENTRALITY_DIVS_SHA256");
    binding.scale_sha256 = env("RJ_AUAU_CENTRALITY_SCALE_SHA256");
    binding.vertex_scale_sha256 =
        env("RJ_AUAU_CENTRALITY_VERTEX_SCALE_SHA256");
    const bool any_path = !binding.divisions.empty() || !binding.scale.empty() ||
                          !binding.vertex_scale.empty();
    const bool any_hash = !binding.divisions_sha256.empty() ||
                          !binding.scale_sha256.empty() ||
                          !binding.vertex_scale_sha256.empty();

    if (!is_auau_data)
    {
      if (!source.empty() || any_path || any_hash)
        detail::bail("RJ_AUAU_CENTRALITY_* must be unset for pp/SIM");
      return binding;
    }
    if (!centrality_will_run)
    {
      if (!source.empty() || any_path || any_hash)
        detail::bail("centrality configured where CentralityReco is skipped");
      return binding;
    }
    if (source != "cdb" && source != "local")
      detail::bail("AuAu DATA requires RJ_AUAU_CENTRALITY_SOURCE=cdb|local");
    if (source == "cdb")
    {
      if (any_path || any_hash)
        detail::bail("CDB centrality mode forbids local override paths/hashes");
      return binding;
    }

    if (run <= 0 || binding.divisions.empty() || binding.scale.empty() ||
        binding.vertex_scale.empty() || binding.divisions_sha256.empty() ||
        binding.scale_sha256.empty() || binding.vertex_scale_sha256.empty())
      detail::bail("local centrality requires run, 3 paths, and 3 expected SHA256s");
    require_run_file(binding.divisions,
                     "cdb_centrality_" + std::to_string(run) + ".root");
    require_run_file(binding.scale,
                     "cdb_centrality_scale_" + std::to_string(run) + ".root");
    require_run_file(binding.vertex_scale,
                     "cdb_centrality_vertex_scale_" + std::to_string(run) + ".root");
    verify_file_sha256(binding.divisions, binding.divisions_sha256);
    verify_file_sha256(binding.scale, binding.scale_sha256);
    verify_file_sha256(binding.vertex_scale, binding.vertex_scale_sha256);
    binding.apply_local = true;
    return binding;
  }

  inline void apply(CentralityReco* centrality, const Binding& binding)
  {
    if (!centrality) detail::bail("null CentralityReco");
    if (binding.apply_local)
    {
      centrality->setOverwriteDivs(binding.divisions);
      centrality->setOverwriteScale(binding.scale);
      centrality->setOverwriteVtx(binding.vertex_scale);
      std::cout << "CENTRALITY_BINDING mode=local domain=AuAu_DATA_MB_ONLY"
                << " run=" << binding.run
                << " divisions=" << binding.divisions
                << " scale=" << binding.scale
                << " vertex_scale=" << binding.vertex_scale
                << " application_count=1"
                << " native_initrun_readback=REQUIRED" << std::endl;
    }
    else
    {
      std::cout << "CENTRALITY_BINDING mode=cdb domain=AuAu_DATA_MB_ONLY"
                << " run=" << binding.run
                << " keys=Centrality,CentralityScale,CentralityVertexScale"
                << " application_count=1" << std::endl;
    }
  }
}

#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
class ScaledTriggerStudyReco : public SubsysReco
{
 public:
  explicit ScaledTriggerStudyReco(const std::string& outFile)
    : SubsysReco("ScaledTriggerStudyReco")
    , m_outFile(outFile)
  {
    if (const char* env = std::getenv("RJ_SCALED_TRIGGER_VZ_MAX_CM"))
    {
      const double v = std::atof(env);
      if (std::isfinite(v) && v > 0.0) m_vzMaxCm = v;
    }
    if (const char* env = std::getenv("RJ_SCALED_TRIGGER_RUNLIST"))
    {
      const std::string s = detail::trim(std::string(env));
      if (!s.empty()) m_runListPath = s;
    }
    if (const char* env = std::getenv("RJ_SCALED_TRIGGER_CENT_STUDY"))
    {
      const std::string s = lower(detail::trim(std::string(env)));
      m_enableCentrality = (s == "1" || s == "true" || s == "yes" || s == "on");
    }
    if (const char* env = std::getenv("RJ_SCALED_TRIGGER_CENT_EDGES"))
    {
      std::vector<int> parsed = parseCentEdges(env);
      if (parsed.size() >= 2)
      {
        m_centEdges = parsed;
        m_enableCentrality = true;
      }
    }
  }

  int Init(PHCompositeNode*) override
  {
    loadRunList();
    if (m_selectedRuns.empty())
    {
      std::cerr << "[ScaledTriggerStudyReco] empty selected-run list: "
                << m_runListPath << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_out = TFile::Open(m_outFile.c_str(), "RECREATE");
    if (!m_out || m_out->IsZombie())
    {
      std::cerr << "[ScaledTriggerStudyReco] could not open output file: "
                << m_outFile << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_hBase = bookHist("MBD_NS_geq_2_vtx_lt_150",
                       "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_MBD_NS_geq_2_vtx_lt_150");
    m_hPho10 = bookHist("Photon_10",
                        "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_10");
    m_hPho12 = bookHist("Photon_12",
                        "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_12");

    if (!m_hBase || !m_hPho10 || !m_hPho12)
    {
      std::cerr << "[ScaledTriggerStudyReco] failed to book required histograms" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    if (m_enableCentrality)
    {
      const size_t nCent = m_centEdges.size() - 1;
      m_hBaseCent.resize(nCent, nullptr);
      m_hPho10Cent.resize(nCent, nullptr);
      m_hPho12Cent.resize(nCent, nullptr);

      for (size_t i = 0; i < nCent; ++i)
      {
        const std::string suffix = centSuffix(i);
        m_hBaseCent[i] = bookHist("MBD_NS_geq_2_vtx_lt_150",
                                  "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_MBD_NS_geq_2_vtx_lt_150" + suffix);
        m_hPho10Cent[i] = bookHist("Photon_10",
                                   "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_10" + suffix);
        m_hPho12Cent[i] = bookHist("Photon_12",
                                   "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_12" + suffix);
        if (!m_hBaseCent[i] || !m_hPho10Cent[i] || !m_hPho12Cent[i])
        {
          std::cerr << "[ScaledTriggerStudyReco] failed to book centrality histograms for "
                    << centLabel(i) << std::endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }
      }
    }

    return Fun4AllReturnCodes::EVENT_OK;
  }

  int process_event(PHCompositeNode* topNode) override
  {
    if (!topNode) return Fun4AllReturnCodes::ABORTEVENT;

    const uint64_t runNumber = currentRun(topNode);
    if (!runNumber) return Fun4AllReturnCodes::EVENT_OK;
    if (!std::binary_search(m_selectedRuns.begin(), m_selectedRuns.end(), runNumber))
      return Fun4AllReturnCodes::EVENT_OK;

    const double vz = recoVz(topNode);
    if (!std::isfinite(vz) || std::fabs(vz) >= m_vzMaxCm)
      return Fun4AllReturnCodes::EVENT_OK;

    Gl1Packet* gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!gl1) gl1 = findNode::getClass<Gl1Packet>(topNode, "14001");
    if (!gl1) return Fun4AllReturnCodes::EVENT_OK;

    const uint64_t scaledVector = gl1->lValue(0, "ScaledVector");
    const float maxEnergy = maxClusterEnergy(topNode);

    const bool passBase = bitIsSet(scaledVector, 14);
    const bool passPho10 = bitIsSet(scaledVector, 22);
    const bool passPho12 = bitIsSet(scaledVector, 23);

    if (passBase) m_hBase->Fill(maxEnergy);
    if (passPho10) m_hPho10->Fill(maxEnergy);
    if (passPho12) m_hPho12->Fill(maxEnergy);

    if (m_enableCentrality)
    {
      const int centIdx = centralityIndex(topNode);
      if (centIdx >= 0)
      {
        if (passBase) m_hBaseCent[centIdx]->Fill(maxEnergy);
        if (passPho10) m_hPho10Cent[centIdx]->Fill(maxEnergy);
        if (passPho12) m_hPho12Cent[centIdx]->Fill(maxEnergy);
      }
    }

    ++m_eventsAccepted;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int End(PHCompositeNode*) override
  {
    if (!m_out || !m_out->IsOpen()) return Fun4AllReturnCodes::EVENT_OK;

    writeHist("MBD_NS_geq_2_vtx_lt_150", m_hBase);
    writeHist("Photon_10", m_hPho10);
    writeHist("Photon_12", m_hPho12);

    if (m_enableCentrality)
    {
      for (size_t i = 0; i < m_hBaseCent.size(); ++i)
      {
        writeHist("MBD_NS_geq_2_vtx_lt_150", m_hBaseCent[i]);
        writeHist("Photon_10", m_hPho10Cent[i]);
        writeHist("Photon_12", m_hPho12Cent[i]);
      }
    }

    m_out->cd();
    TNamed mode("scaledTriggerStudyOnly", "1");
    mode.Write("scaledTriggerStudyOnly", TObject::kOverwrite);
    TNamed centMode("scaledTriggerCentStudy", m_enableCentrality ? "1" : "0");
    centMode.Write("scaledTriggerCentStudy", TObject::kOverwrite);
    TNamed runList("scaledTriggerStudyRunList", m_runListPath.c_str());
    runList.Write("scaledTriggerStudyRunList", TObject::kOverwrite);
    TNamed vzCut("scaledTriggerStudyVzMaxCm", std::to_string(m_vzMaxCm).c_str());
    vzCut.Write("scaledTriggerStudyVzMaxCm", TObject::kOverwrite);
    TNamed centEdges("scaledTriggerStudyCentEdges", joinInts(m_centEdges).c_str());
    centEdges.Write("scaledTriggerStudyCentEdges", TObject::kOverwrite);
    TNamed centMissing("scaledTriggerStudyCentralityMissing",
                       std::to_string(m_centralityMissing).c_str());
    centMissing.Write("scaledTriggerStudyCentralityMissing", TObject::kOverwrite);
    TNamed centInvalid("scaledTriggerStudyCentralityInvalid",
                       std::to_string(m_centralityInvalid).c_str());
    centInvalid.Write("scaledTriggerStudyCentralityInvalid", TObject::kOverwrite);
    TNamed centOut("scaledTriggerStudyCentralityOutOfRange",
                   std::to_string(m_centralityOutOfRange).c_str());
    centOut.Write("scaledTriggerStudyCentralityOutOfRange", TObject::kOverwrite);
    TNamed accepted("scaledTriggerStudyEventsAccepted", std::to_string(m_eventsAccepted).c_str());
    accepted.Write("scaledTriggerStudyEventsAccepted", TObject::kOverwrite);

    m_out->Write("", TObject::kOverwrite);
    m_out->Close();
    delete m_out;
    m_out = nullptr;
    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  static bool bitIsSet(uint64_t vec, int bit)
  {
    return bit >= 0 && bit < 64 && (vec & (1ULL << bit));
  }

  static std::string lower(std::string s)
  {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
  }

  static std::vector<int> parseCentEdges(std::string s)
  {
    for (char& c : s)
    {
      if (c == ',' || c == ';' || c == ':') c = ' ';
    }
    std::istringstream iss(s);
    std::vector<int> edges;
    int edge = 0;
    while (iss >> edge) edges.push_back(edge);
    if (edges.size() < 2) return {};
    if (!std::is_sorted(edges.begin(), edges.end())) return {};
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges.size() >= 2 ? edges : std::vector<int>{};
  }

  static std::string joinInts(const std::vector<int>& values)
  {
    std::ostringstream os;
    for (size_t i = 0; i < values.size(); ++i)
    {
      if (i) os << ",";
      os << values[i];
    }
    return os.str();
  }

  void loadRunList()
  {
    std::ifstream in(m_runListPath);
    uint64_t run = 0;
    while (in >> run)
    {
      if (run > 0) m_selectedRuns.push_back(run);
    }
    std::sort(m_selectedRuns.begin(), m_selectedRuns.end());
    m_selectedRuns.erase(std::unique(m_selectedRuns.begin(), m_selectedRuns.end()),
                         m_selectedRuns.end());
  }

  TH1F* bookHist(const std::string& dirName, const std::string& histName)
  {
    if (!m_out) return nullptr;
    TDirectory* dir = m_out->GetDirectory(dirName.c_str());
    if (!dir) dir = m_out->mkdir(dirName.c_str());
    if (!dir) return nullptr;
    dir->cd();
    TH1F* h = new TH1F(histName.c_str(), "Max Cluster Energy; Cluster Energy [GeV]", 40, 0, 20);
    h->SetDirectory(dir);
    m_out->cd();
    return h;
  }

  std::string centLabel(size_t i) const
  {
    std::ostringstream os;
    os << m_centEdges[i] << "-" << m_centEdges[i + 1] << "%";
    return os.str();
  }

  std::string centSuffix(size_t i) const
  {
    std::ostringstream os;
    os << "_cent" << m_centEdges[i] << "_" << m_centEdges[i + 1];
    return os.str();
  }

  int centralityIndex(PHCompositeNode* topNode)
  {
    CentralityInfo* central = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
    if (!central)
    {
      ++m_centralityMissing;
      return -1;
    }

    const float centile = central->get_centrality_bin(CentralityInfo::PROP::mbd_NS);
    if (!std::isfinite(centile) || centile < 0.0f)
    {
      ++m_centralityInvalid;
      return -1;
    }

    for (size_t i = 0; i + 1 < m_centEdges.size(); ++i)
    {
      if (centile >= static_cast<float>(m_centEdges[i]) &&
          centile < static_cast<float>(m_centEdges[i + 1]))
      {
        return static_cast<int>(i);
      }
    }

    ++m_centralityOutOfRange;
    return -1;
  }

  void writeHist(const std::string& dirName, TH1* h)
  {
    if (!m_out || !h) return;
    TDirectory* dir = m_out->GetDirectory(dirName.c_str());
    if (!dir) dir = m_out->mkdir(dirName.c_str());
    if (!dir) return;
    dir->cd();
    h->Write("", TObject::kOverwrite);
    m_out->cd();
  }

  static uint64_t currentRun(PHCompositeNode* topNode)
  {
    if (auto* evt = findNode::getClass<EventHeader>(topNode, "EventHeader"))
    {
      const int run = evt->get_RunNumber();
      if (run > 0) return static_cast<uint64_t>(run);
    }
    return recoConsts::instance()->get_uint64Flag("TIMESTAMP", 0);
  }

  static double recoVz(PHCompositeNode* topNode)
  {
    if (auto* mbdMap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap"))
    {
      if (!mbdMap->empty())
      {
        MbdVertex* v = mbdMap->begin()->second;
        if (v && std::isfinite(v->get_z())) return v->get_z();
      }
    }
    if (auto* gvMap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap"))
    {
      if (!gvMap->empty())
      {
        GlobalVertex* v = gvMap->begin()->second;
        if (v && std::isfinite(v->get_z())) return v->get_z();
      }
    }
    return std::numeric_limits<double>::quiet_NaN();
  }

  static float maxClusterEnergy(PHCompositeNode* topNode)
  {
    float maxEnergy = 0.0f;
    RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
    if (!clusters) return maxEnergy;

    const auto range = clusters->getClusters();
    for (auto it = range.first; it != range.second; ++it)
    {
      const RawCluster* cl = it->second;
      if (!cl) continue;
      const float e = cl->get_energy();
      if (!std::isfinite(e) || e < 1.0f) continue;
      if (e > maxEnergy) maxEnergy = e;
    }
    return maxEnergy;
  }

  std::string m_outFile;
  // Optional scaled-trigger QA run list; set RJ_SCALED_TRIGGER_RUNLIST to use it.
  std::string m_runListPath;
  std::vector<uint64_t> m_selectedRuns;
  double m_vzMaxCm = 30.0;
  bool m_enableCentrality = false;
  std::vector<int> m_centEdges = {0, 20, 50, 80};
  uint64_t m_centralityMissing = 0;
  uint64_t m_centralityInvalid = 0;
  uint64_t m_centralityOutOfRange = 0;
  uint64_t m_eventsAccepted = 0;
  TFile* m_out = nullptr;
  TH1F* m_hBase = nullptr;
  TH1F* m_hPho10 = nullptr;
  TH1F* m_hPho12 = nullptr;
  std::vector<TH1F*> m_hBaseCent;
  std::vector<TH1F*> m_hPho10Cent;
  std::vector<TH1F*> m_hPho12Cent;
};
#endif

namespace yamlcfg
{
  // Forward declarations (LoadJetRKeys appears before helper definitions below)
  inline std::string ResolveYAMLPath();
  inline bool ReadWholeFile(const std::string& path, std::string& out);
  inline bool StartsWithKey(const std::string& line, const std::string& key);
  inline std::string AfterColon(const std::string& line);
  inline std::string StripQuotes(const std::string& value);
  inline void ParseInlineListDoubles(std::string s, std::vector<double>& out);

  inline std::vector<std::string> LoadJetRKeys(int vlevel)
  {
    std::vector<std::string> outKeys;

    const std::string yamlPath = ResolveYAMLPath();
    std::string yamlText;
    const bool ok = ReadWholeFile(yamlPath, yamlText);

    // Default behavior (baseline): r02 + r04
    std::vector<double> radii = {0.2, 0.4};

    if (ok)
    {
      std::istringstream iss(yamlText);
      for (std::string line; std::getline(iss, line); )
      {
        line = detail::trim(line);
        if (line.empty()) continue;
        if (!line.empty() && line[0] == '#') continue;

        if (StartsWithKey(line, "jet_radii"))
        {
          const std::string rhs = AfterColon(line);
          ParseInlineListDoubles(rhs, radii);
          break;
        }
      }
    }

    auto push_unique = [&](const std::string& k)
    {
      if (std::find(outKeys.begin(), outKeys.end(), k) == outKeys.end())
        outKeys.push_back(k);
    };

    for (double R : radii)
    {
      if (!std::isfinite(R) || R <= 0.0) continue;

      // Map R to "r%02d" where D ~ round(10*R)
      const int D = static_cast<int>(R * 10.0 + 0.5);
      if (D <= 0) continue;

      std::ostringstream rk;
      rk << "r" << std::setw(2) << std::setfill('0') << D;
      push_unique(rk.str());
    }

    if (outKeys.empty())
    {
      outKeys = {"r02","r04"};
      if (vlevel > 0)
        std::cout << "[CFG] jet_radii: no valid entries found -> defaulting to [r02,r04]\n";
    }

    if (vlevel > 0)
    {
      std::cout << "[CFG] jet_radii from YAML (" << yamlPath << "): [";
      for (std::size_t i = 0; i < outKeys.size(); ++i)
        std::cout << outKeys[i] << (i + 1 < outKeys.size() ? ", " : "");
      std::cout << "]\n";
    }

    return outKeys;
  }

  inline bool WantRKey(const std::vector<std::string>& keys, const std::string& key)
  {
    if (keys.empty()) return true;
    return (std::find(keys.begin(), keys.end(), key) != keys.end());
  }
}

namespace yamlcfg
{
    struct Config
    {
        // file provenance
        std::string yamlPath = "";
        std::string yamlText = "";
        
        // baseline defaults (match current behavior)
        double photon_eta_abs_max = 0.7;
        double jet_pt_min = 10.0;
        double back_to_back_dphi_min_pi_fraction = 0.875;
        
        bool   use_vz_cut = true;
        double vz_cut_cm  = 30.0;
        bool   setMinBiasClassifer = false;
        
        std::vector<int> centrality_edges = {0, 10, 20, 40, 60, 80, 100};
        
        bool vertex_reweight_on_pp = true;
        std::string vertex_reweight_file_pp;    // set in analysis_config.yaml
        std::string vertex_reweight_hist_pp = "h_w_iterative";

        bool vertex_reweight_on_auau = false;
        std::string vertex_reweight_file_auau;  // set in analysis_config.yaml
        std::string vertex_reweight_hist_auau = "data_over_MC_ratios/h_zvtx_ratio_data_over_photonJet";
        
        bool centrality_reweight_on = false;
        std::string centrality_reweight_file;   // set in analysis_config.yaml
        std::string centrality_reweight_hist = "nom_cent_rw_hist";
        
        double isoA = 0.490;
        double isoB = 0.037;
        double isoGap = 0.8;
        double isoFixed = 2.0;
        double truthIsoGeV = 4.0;
        double isoConeR = 0.30;
        double isoTowMin = 0.0;
        bool   isSlidingIso = true;
        bool   isSlidingAndFixed = false;
        
        // Per-centrality AuAu sliding WPs: each entry = {aGeV, bPerGeV, sideGapGeV}
        struct CentIsoWP { double aGeV; double bPerGeV; double sideGapGeV; };
        struct PPIsoWP
        {
            double aGeV = 0.0;
            double bPerGeV = 0.0;
            double sideGapGeV = 1.0;
            bool configured = false;
        };
        PPIsoWP ppIsoWPR30;
        PPIsoWP ppIsoWPR40;
        std::vector<CentIsoWP> auauCentIsoWP;
        std::vector<CentIsoWP> auauCentIsoWPR30;
        std::vector<CentIsoWP> auauCentIsoWPR40;
        
        // Photon ID cuts (PPG12 Table 4) baseline
        double pre_e11e33_max = 0.98;
        double pre_et1_min    = 0.60;
        double pre_et1_max    = 1.00;
        double pre_e32e35_min = 0.80;
        double pre_e32e35_max = 1.00;
        double pre_weta_max   = 0.60;
        
        double tight_w_lo            = 0.0;
        double tight_w_hi_intercept  = 0.15;
        double tight_w_hi_slope      = 0.006;
        
        double tight_e11e33_min = 0.40;
        double tight_e11e33_max = 0.98;
        
        double tight_et1_min    = 0.90;
        double tight_et1_max    = 1.00;
        
        double tight_e32e35_min = 0.92;
        double tight_e32e35_max = 1.00;
        
        double pho_dr_max = 0.05;
        double jet_dr_max = 0.3;
        
        std::vector<double> jes3_photon_pt_bins = {15,17,19,21,23,26,35};
        std::vector<double> unfold_reco_photon_pt_bins  = {10,15,17,19,21,23,26,35,40};
        std::vector<double> unfold_truth_photon_pt_bins = {5,10,15,17,19,21,23,26,35,40};
        
        double unfold_jet_pt_start = 0.0;
        double unfold_jet_pt_stop  = 60.0;
        double unfold_jet_pt_step  = 0.5;
        
        std::vector<double> unfold_xj_bins = {0.0,0.20,0.24,0.29,0.35,0.41,0.50,0.60,0.72,0.86,1.03,1.24,1.49,1.78,2.14,3.0};
        std::string leading_response_family = "";
        bool require_towerinfo_truth_matching = false;
        
        // EventDisplay diagnostics payload (EventDisplayTree)
        bool event_display_tree = true;
        int  event_display_tree_max_per_bin = 0;
        std::string clusterUEpipeline = "noSub";

        std::string preselection = "reference";
        std::string tight = "reference";
        std::string nonTight = "reference";

        std::string npb_model_file = "";
        double npb_cut = 0.5;
        std::vector<std::string> npb_features;

        std::string tight_bdt_model_file = "";
        std::string ppg12_base_e_model_file = "";
        double tight_bdt_min_intercept = 0.815625;
        double tight_bdt_min_slope = -0.0015625;
        double tight_bdt_max = 1.0;
        double nontight_bdt_min_intercept = 0.7333333333333333;
        double nontight_bdt_min_slope = -0.01333333333333333;
        double nontight_bdt_max_intercept = 0.684375;
        double nontight_bdt_max_slope = 0.0015625;
        std::vector<std::string> tight_bdt_features;

        std::string auau_npb_model_file = "";
        double auau_npb_cut = 0.5;
        std::vector<std::string> auau_npb_features;

        std::string auau_tight_bdt_model_file = "";
        std::string auau_tight_bdt_expanded_model_dir = "";
        std::string auau_tight_bdt_centINDcontrol_model_file = "";
        std::string auau_tight_bdt_centAsFeat_model_file = "";
        std::vector<std::string> auau_tight_bdt_centDep_model_files;
        std::string auau_tight_bdt_noCent_model_file = "";
        std::string auau_tight_bdt_centInput_model_file = "";
        std::string auau_tight_bdt_centInput3x3_model_file = "";
        std::string auau_tight_bdt_centInputBase3x3_model_file = "";
        std::string auau_tight_bdt_centInputMinOpt_model_file = "";
        std::vector<std::string> auau_tight_bdt_cent3_model_files;
        std::vector<std::string> auau_tight_bdt_cent7_model_files;
        std::vector<std::string> auau_tight_bdt_ptBinCentInput_model_files;
        std::vector<std::string> auau_tight_bdt_ptCent3_model_files;
        std::vector<std::string> auau_tight_bdt_ptCent7_model_files;
        std::vector<std::string> auau_tight_bdt_etFineCentInput_model_files;
        std::vector<std::string> auau_tight_bdt_etFineCent3_model_files;
        std::vector<std::string> auau_tight_bdt_etFineCent7_model_files;
        std::string auau_tight_bdt_etFineCent3_product = "";
        std::string auau_tight_bdt_etFineCent7_product = "";
        std::string auau_tight_bdt_ptBinCentInput_fallback_model_file = "";
        std::vector<std::string> auau_tight_bdt_ptCent3_fallback_model_files;
        std::vector<std::string> auau_tight_bdt_ptCent7_fallback_model_files;
        std::vector<double> auau_tight_bdt_pt_bin_edges;
        std::vector<double> auau_tight_bdt_etfine_pt_bin_edges;
        std::vector<int> auau_tight_bdt_cent3_edges;
        std::vector<int> auau_tight_bdt_cent7_edges;
        double auau_tight_bdt_pt_fallback_min = 35.0;
        double auau_tight_bdt_pt_fallback_max = 40.0;
        double auau_tight_bdt_apply_pt_min = std::numeric_limits<double>::quiet_NaN();
        double auau_tight_bdt_apply_pt_max = std::numeric_limits<double>::quiet_NaN();
        double auau_tight_bdt_min_intercept = 0.8333333333333334;
        double auau_tight_bdt_min_slope = -0.003333333333333336;
        double auau_tight_bdt_max = 1.0;
        std::vector<std::string> auau_tight_bdt_working_point_entries;
        double auau_nontight_bdt_min_intercept = 0.7333333333333333;
        double auau_nontight_bdt_min_slope = -0.01333333333333333;
        double auau_nontight_bdt_max_intercept = 0.6666666666666666;
        double auau_nontight_bdt_max_slope = 0.003333333333333336;
        std::string auau_nontight_bdt_sideband_mode = "etLinear";
        double auau_nontight_bdt_relative_min_offset = -0.20;
        double auau_nontight_bdt_relative_max_offset = -0.03;
        std::vector<std::string> auau_tight_bdt_features;
        std::vector<std::string> auau_tight_bdt_centINDcontrol_features;
        std::vector<std::string> auau_tight_bdt_centAsFeat_features;
        std::vector<std::string> auau_tight_bdt_centAsFeat3x3_features;
        std::vector<std::string> auau_tight_bdt_centAsFeatBase3x3_features;
        std::vector<std::string> auau_tight_bdt_centDep_features;
        std::string auau_tight_mlp_model_dir = "";
        std::string auau_tight_mlp_model_file = "";
        std::string auau_tight_mlp_centInput_model_file = "";
        std::string auau_tight_mlp_noCentBase3x3_model_file = "";
        std::string auau_tight_mlp_centInputBase3x3_model_file = "";
        double auau_tight_mlp_min_intercept = 0.80;
        double auau_tight_mlp_min_slope = 0.0;
        double auau_tight_mlp_max = 1.0;
        double auau_tight_mlp_apply_pt_min = std::numeric_limits<double>::quiet_NaN();
        double auau_tight_mlp_apply_pt_max = std::numeric_limits<double>::quiet_NaN();
        double auau_nontight_mlp_min_intercept = 0.20;
        double auau_nontight_mlp_min_slope = 0.0;
        double auau_nontight_mlp_max_intercept = 0.80;
        double auau_nontight_mlp_max_slope = 0.0;
        std::vector<std::string> auau_tight_mlp_working_point_entries;
        std::string auau_tight_bdt_mlp_stack_model_file = "";
        std::string auau_tight_bdt_mlp_stack_bdt_mode = "auauEtFineCent7BDT";
        std::string auau_tight_bdt_mlp_stack_mlp_mode = "auauCentInputBase3x3MLP";
        double auau_tight_bdt_mlp_stack_min_intercept = 0.80;
        double auau_tight_bdt_mlp_stack_min_slope = 0.0;
        double auau_tight_bdt_mlp_stack_max = 1.0;
        double auau_tight_bdt_mlp_stack_apply_pt_min = std::numeric_limits<double>::quiet_NaN();
        double auau_tight_bdt_mlp_stack_apply_pt_max = std::numeric_limits<double>::quiet_NaN();
        double auau_nontight_bdt_mlp_stack_min_intercept = 0.20;
        double auau_nontight_bdt_mlp_stack_min_slope = 0.0;
        double auau_nontight_bdt_mlp_stack_max_intercept = 0.80;
        double auau_nontight_bdt_mlp_stack_max_slope = 0.0;
        std::vector<std::string> auau_tight_bdt_mlp_stack_working_point_entries;
        std::string auau_tight_logreg_model_file = "";
        double auau_tight_logreg_min_intercept = 0.80;
        double auau_tight_logreg_min_slope = 0.0;
        double auau_tight_logreg_max = 1.0;
        double auau_tight_logreg_apply_pt_min = std::numeric_limits<double>::quiet_NaN();
        double auau_tight_logreg_apply_pt_max = std::numeric_limits<double>::quiet_NaN();
        double auau_nontight_logreg_min_intercept = 0.20;
        double auau_nontight_logreg_min_slope = 0.0;
        double auau_nontight_logreg_max_intercept = 0.80;
        double auau_nontight_logreg_max_slope = 0.0;
        std::vector<std::string> auau_tight_logreg_working_point_entries;

        bool auau_bdt_training_tree = false;
        long long auau_bdt_training_tree_max_entries = 0;
        bool auau_bdt_npb_data_tagging = false;
        double auau_npb_tag_delta_t_cut = -7.0;
        double auau_npb_tag_weta_min = 0.0;
        double auau_npb_tag_away_jet_pt_min = 5.0;
        double auau_npb_tag_away_jet_dphi_min = 1.5707963267948966;
        double auau_npb_tag_time_sample_ns = 17.6;
        double auau_npb_mbd_t0_offset = 0.0;

        bool pp_photonid_extract_only = false;
        bool pp_photonid_training_tree = false;
        bool pp_photonid_ppg12_filter = true;
        bool pp_photonid_require_preselection = false;
        long long pp_photonid_training_tree_max_entries = 0;
        std::string pp_photonid_source_role = "auto";

        bool jet_ml_training_tree = false;
        long long jet_ml_training_tree_max_entries = 0;
        bool jet_ml_correction_enabled = false;
        bool jet_ml_use_as_nominal = false;
        std::string jet_ml_model_file = "";
        std::vector<std::string> jet_ml_features;

        bool doPi0Analysis = false;
    };

    static std::string NormalizePreselectionMode(std::string mode)
    {
        mode = detail::trim(mode);
        std::string key = mode;
        std::transform(key.begin(), key.end(), key.begin(),
                       [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
        if (key == "" || key == "reference") return "reference";
        if (key == "varianta" || key == "newppg12") return "newPPG12";
        if (key == "variantb" || key == "noprecriteria") return "noPreCriteria";
        if (key == "variantc" || key == "onlynpb") return "onlyNPB";
        if (key == "variantd" || key == "refplusnpb") return "refPlusNPB";
        if (key == "variante" || key == "auauonlynpb") return "auauOnlyNPB";
        return mode;
    }

    static bool IsPreselectionMode(const std::string& mode)
    {
        return mode == "reference" || mode == "newPPG12" || mode == "noPreCriteria" ||
               mode == "onlyNPB" || mode == "refPlusNPB" || mode == "auauOnlyNPB";
    }

    static bool PreselectionUsesNPB(const std::string& mode)
    {
        return mode == "newPPG12" || mode == "onlyNPB" || mode == "refPlusNPB";
    }

    static bool PreselectionUsesAuAuNPB(const std::string& mode)
    {
        return mode == "auauOnlyNPB";
    }

    static std::string NormalizeTightMode(std::string mode)
    {
        mode = detail::trim(mode);
        std::string key = mode;
        std::transform(key.begin(), key.end(), key.begin(),
                       [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
        if (key == "" || key == "reference") return "reference";
        if (key == "varianta" || key == "newppg12") return "newPPG12";
        if (key == "variantb" || key == "auauembeddedbdt") return "auauEmbeddedBDT";
        if (key == "centindcontrol") return "centINDcontrol";
        if (key == "centasfeat" || key == "centasfeature") return "centAsFeat";
        if (key == "centdepbdts" || key == "centdepbdt") return "centDepBDTs";
        if (key == "auaunocentbdt") return "auauNoCentBDT";
        if (key == "auaucentinputbdt") return "auauCentInputBDT";
        if (key == "auaucentinput3x3bdt") return "auauCentInput3x3BDT";
        if (key == "auaucentinputbase3x3bdt") return "auauCentInputBase3x3BDT";
        if (key == "auaucentinputminoptbdt") return "auauCentInputMinOptBDT";
        if (key == "auaucent3bdt") return "auauCent3BDT";
        if (key == "auaucent7bdt") return "auauCent7BDT";
        if (key == "auauptbincentinputbdt") return "auauPtBinCentInputBDT";
        if (key == "auauptcent3bdt") return "auauPtCent3BDT";
        if (key == "auauptcent7bdt") return "auauPtCent7BDT";
        if (key == "auauetfinecentinputbdt") return "auauEtFineCentInputBDT";
        if (key == "auauetfinecent3bdt") return "auauEtFineCent3BDT";
        if (key == "auauetfinecent7bdt") return "auauEtFineCent7BDT";
        if (key == "auaucentinputmlp") return "auauCentInputMLP";
        if (key == "auaunocentbase3x3mlp") return "auauNoCentBase3x3MLP";
        if (key == "auaucentinputbase3x3mlp") return "auauCentInputBase3x3MLP";
        if (key == "auauhighptdistilledkitchenmlp") return "auauHighPtDistilledKitchenMLP";
        if (key == "auaubdtmlpstack" || key == "auaustackbdtmlp" || key == "bdtmlpstack") return "auauBDTMLPStack";
        if (key == "auautightlogreg" || key == "auaulogreg" || key == "logreg") return "auauTightLogReg";
        return mode;
    }

    static std::string NormalizeNonTightMode(std::string mode)
    {
        mode = detail::trim(mode);
        std::string key = mode;
        std::transform(key.begin(), key.end(), key.begin(),
                       [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
        if (key == "" || key == "reference") return "reference";
        if (key == "varianta" || key == "newppg12" || key == "bdtsideband") return "newPPG12";
        if (key == "variantb" || key == "auaubdtsideband") return "auauBDTSideband";
        if (key == "variantc" || key == "auaubdtcomplement") return "auauBDTComplement";
        if (key == "auaumlpsideband") return "auauMLPSideband";
        if (key == "auaumlpcomplement") return "auauMLPComplement";
        if (key == "auaubdtmlpstacksideband") return "auauBDTMLPStackSideband";
        if (key == "auaubdtmlpstackcomplement") return "auauBDTMLPStackComplement";
        if (key == "auaulogregsideband") return "auauLogRegSideband";
        if (key == "auaulogregcomplement") return "auauLogRegComplement";
        if (key == "centindcontrol") return "centINDcontrol";
        if (key == "centasfeat" || key == "centasfeature") return "centAsFeat";
        if (key == "centdepbdts" || key == "centdepbdt") return "centDepBDTs";
        return mode;
    }

    static bool IsTightMode(const std::string& mode)
    {
        return mode == "reference" || mode == "newPPG12" || mode == "auauEmbeddedBDT" ||
               mode == "centINDcontrol" || mode == "centAsFeat" || mode == "centDepBDTs" ||
               mode == "auauNoCentBDT" || mode == "auauCentInputBDT" ||
               mode == "auauCentInput3x3BDT" || mode == "auauCentInputBase3x3BDT" ||
               mode == "auauCentInputMinOptBDT" || mode == "auauCent3BDT" ||
               mode == "auauCent7BDT" || mode == "auauPtBinCentInputBDT" ||
               mode == "auauPtCent3BDT" || mode == "auauPtCent7BDT" ||
               mode == "auauEtFineCentInputBDT" || mode == "auauEtFineCent3BDT" ||
               mode == "auauEtFineCent7BDT" || mode == "auauCentInputMLP" ||
               mode == "auauNoCentBase3x3MLP" || mode == "auauCentInputBase3x3MLP" ||
               mode == "auauHighPtDistilledKitchenMLP" || mode == "auauBDTMLPStack" ||
               mode == "auauTightLogReg";
    }

    static bool IsNonTightMode(const std::string& mode)
    {
        return mode == "reference" || mode == "newPPG12" || mode == "auauBDTSideband" ||
               mode == "auauBDTComplement" || mode == "auauMLPSideband" ||
               mode == "auauMLPComplement" || mode == "auauBDTMLPStackSideband" ||
               mode == "auauBDTMLPStackComplement" || mode == "auauLogRegSideband" ||
               mode == "auauLogRegComplement" || mode == "centINDcontrol" ||
               mode == "centAsFeat" || mode == "centDepBDTs";
    }

    static bool IsAuAuTightBDTMode(const std::string& mode)
    {
        return mode == "auauEmbeddedBDT" || mode == "centINDcontrol" ||
               mode == "centAsFeat" || mode == "centDepBDTs" ||
               mode == "auauNoCentBDT" || mode == "auauCentInputBDT" ||
               mode == "auauCentInput3x3BDT" || mode == "auauCentInputBase3x3BDT" ||
               mode == "auauCentInputMinOptBDT" || mode == "auauCent3BDT" ||
               mode == "auauCent7BDT" || mode == "auauPtBinCentInputBDT" ||
               mode == "auauPtCent3BDT" || mode == "auauPtCent7BDT" ||
               mode == "auauEtFineCentInputBDT" || mode == "auauEtFineCent3BDT" ||
               mode == "auauEtFineCent7BDT";
    }

    static bool IsAuAuTightMLPMode(const std::string& mode)
    {
        return mode == "auauCentInputMLP" ||
               mode == "auauNoCentBase3x3MLP" ||
               mode == "auauCentInputBase3x3MLP" ||
               mode == "auauHighPtDistilledKitchenMLP";
    }

    static bool IsAuAuTightBDTMLPStackMode(const std::string& mode)
    {
        return mode == "auauBDTMLPStack";
    }

    static bool IsAuAuTightLogRegMode(const std::string& mode)
    {
        return mode == "auauTightLogReg";
    }

    inline std::string DefaultYAMLPath()
    {
        std::string here = __FILE__;
        const std::size_t p = here.find_last_of('/');
        if (p == std::string::npos) return std::string("analysis_config.yaml");
        return here.substr(0, p) + "/analysis_config.yaml";
    }
    
    inline std::string ResolveYAMLPath()
    {
        if (const char* env = std::getenv("RJ_CONFIG_YAML"))
        {
            std::string s = detail::trim(std::string(env));
            if (!s.empty()) return s;
        }
        return DefaultYAMLPath();
    }
    
    inline bool ReadWholeFile(const std::string& path, std::string& out)
    {
        std::ifstream in(path);
        if (!in.is_open()) return false;
        std::ostringstream ss;
        ss << in.rdbuf();
        out = ss.str();
        return true;
    }
    
    inline bool StartsWithKey(const std::string& line, const std::string& key)
    {
        const std::string k = key + ":";
        if (line.size() < k.size()) return false;
        return (line.compare(0, k.size(), k) == 0);
    }
    
    inline std::string AfterColon(const std::string& line)
    {
        const std::size_t c = line.find(':');
        if (c == std::string::npos) return std::string{};
        return detail::trim(line.substr(c + 1));
    }

    inline std::string StripQuotes(const std::string& value)
    {
        std::string out = detail::trim(value);
        if (out.size() >= 2)
        {
            const char first = out.front();
            const char last = out.back();
            if ((first == '"' && last == '"') || (first == '\'' && last == '\''))
            {
                out = out.substr(1, out.size() - 2);
            }
        }
        return out;
    }
    
    inline bool ParseBool(std::string s, bool& out)
    {
        s = detail::trim(s);
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
        if (s == "true"  || s == "1" || s == "yes" || s == "on")  { out = true;  return true; }
        if (s == "false" || s == "0" || s == "no"  || s == "off") { out = false; return true; }
        return false;
    }
    
    inline bool ParseDouble(const std::string& s, double& out)
    {
        try { out = std::stod(detail::trim(s)); return true; } catch (...) { return false; }
    }
    
    inline bool ParseInt(const std::string& s, int& out)
    {
        try { out = std::stoi(detail::trim(s)); return true; } catch (...) { return false; }
    }
    
    inline void ParseInlineListDoubles(std::string s, std::vector<double>& out)
    {
        out.clear();
        s = detail::trim(s);
        const std::size_t l = s.find('[');
        const std::size_t r = s.find(']');
        if (l == std::string::npos || r == std::string::npos || r <= l) return;
        std::string inner = s.substr(l + 1, r - l - 1);
        std::stringstream ss(inner);
        std::string tok;
        while (std::getline(ss, tok, ','))
        {
            double v = 0.0;
            if (ParseDouble(tok, v)) out.push_back(v);
        }
    }
    
    inline void ParseInlineListInts(std::string s, std::vector<int>& out)
    {
        out.clear();
        s = detail::trim(s);
        const std::size_t l = s.find('[');
        const std::size_t r = s.find(']');
        if (l == std::string::npos || r == std::string::npos || r <= l) return;
        std::string inner = s.substr(l + 1, r - l - 1);
        std::stringstream ss(inner);
        std::string tok;
        while (std::getline(ss, tok, ','))
        {
            int v = 0;
            if (ParseInt(tok, v)) out.push_back(v);
        }
    }

    inline void ParseInlineListStrings(std::string s, std::vector<std::string>& out)
    {
        out.clear();
        s = detail::trim(s);
        const std::size_t l = s.find('[');
        const std::size_t r = s.find(']');
        if (l == std::string::npos || r == std::string::npos || r <= l) return;
        std::string inner = s.substr(l + 1, r - l - 1);
        std::stringstream ss(inner);
        std::string tok;
        while (std::getline(ss, tok, ','))
        {
            tok = detail::trim(tok);
            if (tok.size() >= 2)
            {
                const char q0 = tok.front();
                const char q1 = tok.back();
                if ((q0 == '"' && q1 == '"') || (q0 == '\'' && q1 == '\''))
                {
                    tok = tok.substr(1, tok.size() - 2);
                }
            }
            if (!tok.empty()) out.push_back(tok);
        }
    }
    
    inline void ParseInlineMapDoubles(std::string s, std::map<std::string, double>& out)
    {
        out.clear();
        s = detail::trim(s);
        const std::size_t l = s.find('{');
        const std::size_t r = s.find('}');
        if (l == std::string::npos || r == std::string::npos || r <= l) return;
        std::string inner = s.substr(l + 1, r - l - 1);
        std::stringstream ss(inner);
        std::string pair;
        while (std::getline(ss, pair, ','))
        {
            const std::size_t c = pair.find(':');
            if (c == std::string::npos) continue;
            std::string k = detail::trim(pair.substr(0, c));
            std::string v = detail::trim(pair.substr(c + 1));
            double dv = 0.0;
            if (ParseDouble(v, dv)) out[k] = dv;
        }
    }
    
    inline Config LoadConfig()
    {
        Config cfg;
        cfg.yamlPath = ResolveYAMLPath();
        
        const bool okRead = ReadWholeFile(cfg.yamlPath, cfg.yamlText);
        
        // Local verbosity for YAML parsing diagnostics (independent of the later banner)
        int yamlV = 0;
        if (const char* venv = std::getenv("RJ_VERBOSITY"))
        {
            try { yamlV = std::stoi(venv); } catch (...) { yamlV = 0; }
        }
        
        bool strict = false;
        if (const char* env = std::getenv("RJ_YAML_STRICT")) strict = (std::atoi(env) != 0);
        
        auto warn_parse = [&](const std::string& key, const std::string& rhs, const std::string& why)
        {
            if (yamlV > 0)
            {
                std::cout << "[CFG][WARN] YAML parse failed for key '" << key << "'"
                << " (rhs='" << rhs << "')"
                << " -> " << why
                << " (keeping default)\n";
            }
            if (strict)
            {
                std::ostringstream oss;
                oss << "YAML parse failed for key '" << key << "' (rhs='" << rhs << "'): " << why;
                throw std::runtime_error(oss.str());
            }
        };
        
        auto info_parse = [&](const std::string& msg)
        {
            if (yamlV > 0) std::cout << msg << "\n";
        };
        
        if (!okRead || cfg.yamlText.empty())
        {
            if (yamlV > 0)
            {
                std::cout << "[CFG][WARN] Could not read YAML (or file empty): " << cfg.yamlPath
                << " (all values will remain defaults)\n";
            }
            if (strict)
            {
                throw std::runtime_error(std::string("Could not read YAML: ") + cfg.yamlPath);
            }
            return cfg;
        }
        
        std::string yamlSection;
        std::istringstream iss(cfg.yamlText);
        std::string pendingLine;
        bool hasPendingLine = false;

        auto parse_cent_iso_wp_list = [&](const std::string& key, std::vector<Config::CentIsoWP>& out)
        {
            out.clear();
            std::string line;
            while (std::getline(iss, line))
            {
                line = detail::trim(line);
                if (line.empty()) continue;
                if (line[0] == '#') continue;
                if (line[0] != '-')
                {
                    pendingLine = line;
                    hasPendingLine = true;
                    break;
                }
                std::map<std::string, double> m;
                ParseInlineMapDoubles(line.substr(1), m);
                Config::CentIsoWP wp{};
                wp.aGeV       = m.count("aGeV")       ? m["aGeV"]       : cfg.isoA;
                wp.bPerGeV    = m.count("bPerGeV")     ? m["bPerGeV"]    : cfg.isoB;
                wp.sideGapGeV = m.count("sideGapGeV")  ? m["sideGapGeV"] : cfg.isoGap;
                out.push_back(wp);
            }
            if (yamlV > 0)
            {
                std::cout << "[CFG] " << key << ": parsed " << out.size() << " entries\n";
            }
        };

        auto yaml_has_key = [&](const std::string& key) -> bool
        {
            std::istringstream text(cfg.yamlText);
            std::string raw;
            while (std::getline(text, raw))
            {
                std::string line = detail::trim(raw);
                if (line.empty() || line[0] == '#') continue;
                if (StartsWithKey(line, key)) return true;
            }
            return false;
        };

        auto recover_cent_iso_wp_block = [&](const std::string& key, std::vector<Config::CentIsoWP>& out)
        {
            out.clear();
            std::istringstream text(cfg.yamlText);
            std::string raw;
            bool inBlock = false;
            while (std::getline(text, raw))
            {
                std::string line = detail::trim(raw);
                if (line.empty() || line[0] == '#') continue;
                if (!inBlock)
                {
                    if (StartsWithKey(line, key)) inBlock = true;
                    continue;
                }
                if (line[0] != '-') break;
                std::map<std::string, double> m;
                ParseInlineMapDoubles(line.substr(1), m);
                Config::CentIsoWP wp{};
                wp.aGeV       = m.count("aGeV")       ? m["aGeV"]       : cfg.isoA;
                wp.bPerGeV    = m.count("bPerGeV")    ? m["bPerGeV"]    : cfg.isoB;
                wp.sideGapGeV = m.count("sideGapGeV") ? m["sideGapGeV"] : cfg.isoGap;
                out.push_back(wp);
            }
            if (yamlV > 0 && !out.empty())
            {
                std::cout << "[CFG] " << key << ": recovered " << out.size()
                          << " entries from YAML block\n";
            }
        };

        for (std::string rawLine; ; )
        {
            if (hasPendingLine)
            {
                rawLine = pendingLine;
                pendingLine.clear();
                hasPendingLine = false;
            }
            else if (!std::getline(iss, rawLine))
            {
                break;
            }
            std::string line = detail::trim(rawLine);
            if (line.empty()) continue;
            if (!line.empty() && line[0] == '#') continue;

            if (line.size() > 1 && line.back() == ':' && line.find('[') == std::string::npos && line.find('{') == std::string::npos)
            {
                yamlSection = detail::trim(line.substr(0, line.size() - 1));
                continue;
            }

            if (StartsWithKey(line, "photon_eta_abs_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.photon_eta_abs_max))
                    warn_parse("photon_eta_abs_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "jet_pt_min"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.jet_pt_min))
                {
                    std::vector<double> v;
                    ParseInlineListDoubles(rhs, v);
                    if (!v.empty())
                    {
                        cfg.jet_pt_min = v.front();
                        std::ostringstream oss;
                        oss << "[CFG] jet_pt_min is a list (n=" << v.size() << "); using first value = " << cfg.jet_pt_min;
                        info_parse(oss.str());
                    }
                    else
                    {
                        warn_parse("jet_pt_min", rhs, "expected a scalar double or an inline list [..]");
                    }
                }
            }
            else if (StartsWithKey(line, "back_to_back_dphi_min_pi_fraction"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.back_to_back_dphi_min_pi_fraction))
                {
                    std::vector<double> v;
                    ParseInlineListDoubles(rhs, v);
                    if (!v.empty())
                    {
                        cfg.back_to_back_dphi_min_pi_fraction = v.front();
                        std::ostringstream oss;
                        oss << "[CFG] back_to_back_dphi_min_pi_fraction is a list (n=" << v.size()
                        << "); using first value = " << cfg.back_to_back_dphi_min_pi_fraction;
                        info_parse(oss.str());
                    }
                    else
                    {
                        warn_parse("back_to_back_dphi_min_pi_fraction", rhs, "expected a scalar double or an inline list [..]");
                    }
                }
            }
            else if (StartsWithKey(line, "use_vz_cut"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.use_vz_cut))
                    warn_parse("use_vz_cut", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "vz_cut_cm"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.vz_cut_cm))
                {
                    // Every other inline-list key already accepts this
                    // form.  A single value is the expanded-capture
                    // spelling; a multi-entry list is genuinely
                    // ambiguous and keeps the default rather than
                    // silently picking one.  The warning is never
                    // verbosity-gated: an unnoticed fallback here
                    // narrows the entire captured vertex range.
                    std::vector<double> parsedVzCut;
                    ParseInlineListDoubles(rhs, parsedVzCut);
                    if (parsedVzCut.size() == 1 &&
                        std::isfinite(parsedVzCut.front()) && parsedVzCut.front() > 0.0)
                    {
                        cfg.vz_cut_cm = parsedVzCut.front();
                    }
                    else
                    {
                        std::cout << "[CFG][WARN] vz_cut_cm could not be parsed (rhs='"
                                  << rhs << "'); keeping " << cfg.vz_cut_cm
                                  << " cm. Give one value, not a multi-entry list."
                                  << std::endl;
                    }
                }
            }
            else if (StartsWithKey(line, "setMinBiasClassifer") || StartsWithKey(line, "setMinBiasClassifier"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.setMinBiasClassifer))
                    warn_parse("setMinBiasClassifer", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "coneR"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.isoConeR))
                {
                    std::vector<double> v;
                    ParseInlineListDoubles(rhs, v);
                    if (!v.empty())
                    {
                        cfg.isoConeR = v.front();
                        std::ostringstream oss;
                        oss << "[CFG] coneR is a list (n=" << v.size() << "); using first value = " << cfg.isoConeR;
                        info_parse(oss.str());
                    }
                    else
                    {
                        warn_parse("coneR", rhs, "expected a scalar double or an inline list [..]");
                    }
                }
            }
            else if (StartsWithKey(line, "centrality_edges"))
            {
              const std::string rhs = AfterColon(line);
              std::vector<int> v;
              ParseInlineListInts(rhs, v);
              if (v.size() >= 2)
              {
                cfg.centrality_edges = v;
              }
              else
              {
                warn_parse("centrality_edges", rhs, "expected an inline list of integers with size >= 2 (e.g. [0, 10, 20, 40, 60, 80, 100])");
              }
            }
            else if (StartsWithKey(line, "vertex_reweight_on_pp"))
            {
              const std::string rhs = AfterColon(line);
              if (!ParseBool(rhs, cfg.vertex_reweight_on_pp))
                warn_parse("vertex_reweight_on_pp", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "vertex_reweight_file_pp"))
            {
              cfg.vertex_reweight_file_pp = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "vertex_reweight_hist_pp"))
            {
              cfg.vertex_reweight_hist_pp = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "vertex_reweight_on_auau"))
            {
              const std::string rhs = AfterColon(line);
              if (!ParseBool(rhs, cfg.vertex_reweight_on_auau))
                warn_parse("vertex_reweight_on_auau", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "vertex_reweight_file_auau"))
            {
              cfg.vertex_reweight_file_auau = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "vertex_reweight_hist_auau"))
            {
              cfg.vertex_reweight_hist_auau = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "vertex_reweight_on"))
            {
              const std::string rhs = AfterColon(line);
              if (!ParseBool(rhs, cfg.vertex_reweight_on_auau))
                warn_parse("vertex_reweight_on", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "vertex_reweight_file"))
            {
              cfg.vertex_reweight_file_auau = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "vertex_reweight_hist"))
            {
              cfg.vertex_reweight_hist_auau = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "centrality_reweight_on"))
            {
              const std::string rhs = AfterColon(line);
              if (!ParseBool(rhs, cfg.centrality_reweight_on))
                warn_parse("centrality_reweight_on", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "centrality_reweight_file"))
            {
              cfg.centrality_reweight_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "centrality_reweight_hist"))
            {
              cfg.centrality_reweight_hist = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "clusterUEpipeline"))
            {
                std::string rhs = AfterColon(line);
                std::vector<std::string> vals;
                ParseInlineListStrings(rhs, vals);
                if (!vals.empty()) rhs = vals.front();
                rhs = detail::trim(rhs);
                if (rhs == "true" || rhs == "1")
                    cfg.clusterUEpipeline = "variantA";
                else if (rhs == "false" || rhs == "0")
                    cfg.clusterUEpipeline = "noSub";
                else if (rhs == "noSub" || rhs == "baseVariant" || rhs == "variantA" || rhs == "variantB")
                    cfg.clusterUEpipeline = rhs;
                else
                    warn_parse("clusterUEpipeline", rhs, "expected noSub|baseVariant|variantA|variantB (or inline list [..])");
            }
            else if (StartsWithKey(line, "preselection"))
            {
                std::string rhs = AfterColon(line);
                std::vector<std::string> vals;
                ParseInlineListStrings(rhs, vals);
                if (!vals.empty()) rhs = vals.front();
                rhs = NormalizePreselectionMode(rhs);
                if (IsPreselectionMode(rhs))
                    cfg.preselection = rhs;
                else
                    warn_parse("preselection", rhs, "expected reference|newPPG12|noPreCriteria|onlyNPB|refPlusNPB|auauOnlyNPB (old variantA/B/C/D/E aliases also accepted)");
            }
            else if (StartsWithKey(line, "tight"))
            {
                std::string rhs = AfterColon(line);
                std::vector<std::string> vals;
                ParseInlineListStrings(rhs, vals);
                if (!vals.empty()) rhs = vals.front();
                rhs = NormalizeTightMode(rhs);
                if (IsTightMode(rhs))
                    cfg.tight = rhs;
                else
                    warn_parse("tight", rhs, "expected reference|newPPG12|auauEmbeddedBDT|centINDcontrol|centAsFeat|centDepBDTs or an auau* validation BDT mode");
            }
            else if (StartsWithKey(line, "nonTight"))
            {
                std::string rhs = AfterColon(line);
                std::vector<std::string> vals;
                ParseInlineListStrings(rhs, vals);
                if (!vals.empty()) rhs = vals.front();
                rhs = NormalizeNonTightMode(rhs);
                if (IsNonTightMode(rhs))
                    cfg.nonTight = rhs;
                else
                    warn_parse("nonTight", rhs, "expected reference|newPPG12|auauBDTSideband|auauBDTComplement|centINDcontrol|centAsFeat|centDepBDTs");
            }
            else if (StartsWithKey(line, "npb_model_file"))
            {
                cfg.npb_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "npb_cut"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.npb_cut))
                    warn_parse("npb_cut", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "npb_features"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.npb_features);
                if (cfg.npb_features.empty())
                    warn_parse("npb_features", rhs, "expected an inline list of feature names");
            }
            else if (StartsWithKey(line, "tight_bdt_model_file"))
            {
                cfg.tight_bdt_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "ppg12_base_e_model_file"))
            {
                cfg.ppg12_base_e_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "tight_bdt_min_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.tight_bdt_min_intercept))
                    warn_parse("tight_bdt_min_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "tight_bdt_min_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.tight_bdt_min_slope))
                    warn_parse("tight_bdt_min_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "tight_bdt_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.tight_bdt_max))
                    warn_parse("tight_bdt_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "nontight_bdt_min_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.nontight_bdt_min_intercept))
                    warn_parse("nontight_bdt_min_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "nontight_bdt_min_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.nontight_bdt_min_slope))
                    warn_parse("nontight_bdt_min_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "nontight_bdt_max_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.nontight_bdt_max_intercept))
                    warn_parse("nontight_bdt_max_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "nontight_bdt_max_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.nontight_bdt_max_slope))
                    warn_parse("nontight_bdt_max_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "tight_bdt_features"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.tight_bdt_features);
                if (cfg.tight_bdt_features.empty())
                    warn_parse("tight_bdt_features", rhs, "expected an inline list of feature names");
            }
            else if (StartsWithKey(line, "auau_npb_model_file"))
            {
                cfg.auau_npb_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_npb_cut"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_npb_cut))
                    warn_parse("auau_npb_cut", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_npb_features"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_npb_features);
                if (cfg.auau_npb_features.empty())
                    warn_parse("auau_npb_features", rhs, "expected an inline list of feature names");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_model_file"))
            {
                cfg.auau_tight_bdt_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_expanded_model_dir"))
            {
                cfg.auau_tight_bdt_expanded_model_dir = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centINDcontrol_model_file"))
            {
                cfg.auau_tight_bdt_centINDcontrol_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centAsFeat_model_file"))
            {
                cfg.auau_tight_bdt_centAsFeat_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centDep_model_files"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_centDep_model_files);
                if (cfg.auau_tight_bdt_centDep_model_files.empty())
                    warn_parse("auau_tight_bdt_centDep_model_files", rhs, "expected an inline list of model ROOT files");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_noCent_model_file"))
            {
                cfg.auau_tight_bdt_noCent_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centInput_model_file"))
            {
                cfg.auau_tight_bdt_centInput_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centInput3x3_model_file"))
            {
                cfg.auau_tight_bdt_centInput3x3_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centInputBase3x3_model_file"))
            {
                cfg.auau_tight_bdt_centInputBase3x3_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centInputMinOpt_model_file"))
            {
                cfg.auau_tight_bdt_centInputMinOpt_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_cent3_model_files"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_cent3_model_files);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_cent7_model_files"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_cent7_model_files);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_ptBinCentInput_model_files"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_ptBinCentInput_model_files);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_ptCent3_model_files"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_ptCent3_model_files);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_ptCent7_model_files"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_ptCent7_model_files);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_etFineCentInput_model_files"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_etFineCentInput_model_files);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_etFineCent3_model_files"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_etFineCent3_model_files);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_etFineCent7_model_files"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_etFineCent7_model_files);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_etFineCent3_product"))
            {
                cfg.auau_tight_bdt_etFineCent3_product = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_etFineCent7_product"))
            {
                cfg.auau_tight_bdt_etFineCent7_product = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_ptBinCentInput_fallback_model_file"))
            {
                cfg.auau_tight_bdt_ptBinCentInput_fallback_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_ptCent3_fallback_model_files"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_ptCent3_fallback_model_files);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_ptCent7_fallback_model_files"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_ptCent7_fallback_model_files);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_pt_bin_edges"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListDoubles(rhs, cfg.auau_tight_bdt_pt_bin_edges);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_etfine_pt_bin_edges"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListDoubles(rhs, cfg.auau_tight_bdt_etfine_pt_bin_edges);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_cent3_edges"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListInts(rhs, cfg.auau_tight_bdt_cent3_edges);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_cent7_edges"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListInts(rhs, cfg.auau_tight_bdt_cent7_edges);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_pt_fallback_min"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_pt_fallback_min))
                    warn_parse("auau_tight_bdt_pt_fallback_min", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_pt_fallback_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_pt_fallback_max))
                    warn_parse("auau_tight_bdt_pt_fallback_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_apply_pt_min"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_apply_pt_min))
                    warn_parse("auau_tight_bdt_apply_pt_min", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_apply_pt_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_apply_pt_max))
                    warn_parse("auau_tight_bdt_apply_pt_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_min_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_min_intercept))
                    warn_parse("auau_tight_bdt_min_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_min_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_min_slope))
                    warn_parse("auau_tight_bdt_min_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_max))
                    warn_parse("auau_tight_bdt_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_working_point_entries"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_working_point_entries);
            }
            else if (StartsWithKey(line, "auau_nontight_bdt_min_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_bdt_min_intercept))
                    warn_parse("auau_nontight_bdt_min_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_bdt_min_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_bdt_min_slope))
                    warn_parse("auau_nontight_bdt_min_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_bdt_max_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_bdt_max_intercept))
                    warn_parse("auau_nontight_bdt_max_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_bdt_max_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_bdt_max_slope))
                    warn_parse("auau_nontight_bdt_max_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_bdt_sideband_mode"))
            {
                cfg.auau_nontight_bdt_sideband_mode = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_nontight_bdt_relative_min_offset"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_bdt_relative_min_offset))
                    warn_parse("auau_nontight_bdt_relative_min_offset", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_bdt_relative_max_offset"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_bdt_relative_max_offset))
                    warn_parse("auau_nontight_bdt_relative_max_offset", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_features"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_features);
                if (cfg.auau_tight_bdt_features.empty())
                    warn_parse("auau_tight_bdt_features", rhs, "expected an inline list of feature names");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centINDcontrol_features"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_centINDcontrol_features);
                if (cfg.auau_tight_bdt_centINDcontrol_features.empty())
                    warn_parse("auau_tight_bdt_centINDcontrol_features", rhs, "expected an inline list of feature names");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centAsFeat_features"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_centAsFeat_features);
                if (cfg.auau_tight_bdt_centAsFeat_features.empty())
                    warn_parse("auau_tight_bdt_centAsFeat_features", rhs, "expected an inline list of feature names");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centAsFeat3x3_features"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_centAsFeat3x3_features);
                if (cfg.auau_tight_bdt_centAsFeat3x3_features.empty())
                    warn_parse("auau_tight_bdt_centAsFeat3x3_features", rhs, "expected an inline list of feature names");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centAsFeatBase3x3_features"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_centAsFeatBase3x3_features);
                if (cfg.auau_tight_bdt_centAsFeatBase3x3_features.empty())
                    warn_parse("auau_tight_bdt_centAsFeatBase3x3_features", rhs, "expected an inline list of feature names");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_centDep_features"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_centDep_features);
                if (cfg.auau_tight_bdt_centDep_features.empty())
                    warn_parse("auau_tight_bdt_centDep_features", rhs, "expected an inline list of feature names");
            }
            else if (StartsWithKey(line, "auau_tight_mlp_model_dir"))
            {
                cfg.auau_tight_mlp_model_dir = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_mlp_model_file"))
            {
                cfg.auau_tight_mlp_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_mlp_centInput_model_file"))
            {
                cfg.auau_tight_mlp_centInput_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_mlp_noCentBase3x3_model_file"))
            {
                cfg.auau_tight_mlp_noCentBase3x3_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_mlp_centInputBase3x3_model_file"))
            {
                cfg.auau_tight_mlp_centInputBase3x3_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_mlp_min_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_mlp_min_intercept))
                    warn_parse("auau_tight_mlp_min_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_mlp_min_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_mlp_min_slope))
                    warn_parse("auau_tight_mlp_min_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_mlp_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_mlp_max))
                    warn_parse("auau_tight_mlp_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_mlp_apply_pt_min"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_mlp_apply_pt_min))
                    warn_parse("auau_tight_mlp_apply_pt_min", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_mlp_apply_pt_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_mlp_apply_pt_max))
                    warn_parse("auau_tight_mlp_apply_pt_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_mlp_min_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_mlp_min_intercept))
                    warn_parse("auau_nontight_mlp_min_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_mlp_min_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_mlp_min_slope))
                    warn_parse("auau_nontight_mlp_min_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_mlp_max_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_mlp_max_intercept))
                    warn_parse("auau_nontight_mlp_max_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_mlp_max_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_mlp_max_slope))
                    warn_parse("auau_nontight_mlp_max_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_mlp_working_point_entries"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_mlp_working_point_entries);
            }
            else if (StartsWithKey(line, "auau_tight_bdt_mlp_stack_model_file"))
            {
                cfg.auau_tight_bdt_mlp_stack_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_mlp_stack_bdt_mode"))
            {
                cfg.auau_tight_bdt_mlp_stack_bdt_mode = NormalizeTightMode(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_mlp_stack_mlp_mode"))
            {
                cfg.auau_tight_bdt_mlp_stack_mlp_mode = NormalizeTightMode(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_bdt_mlp_stack_min_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_mlp_stack_min_intercept))
                    warn_parse("auau_tight_bdt_mlp_stack_min_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_mlp_stack_min_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_mlp_stack_min_slope))
                    warn_parse("auau_tight_bdt_mlp_stack_min_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_mlp_stack_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_mlp_stack_max))
                    warn_parse("auau_tight_bdt_mlp_stack_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_mlp_stack_apply_pt_min"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_mlp_stack_apply_pt_min))
                    warn_parse("auau_tight_bdt_mlp_stack_apply_pt_min", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_mlp_stack_apply_pt_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_bdt_mlp_stack_apply_pt_max))
                    warn_parse("auau_tight_bdt_mlp_stack_apply_pt_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_bdt_mlp_stack_min_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_bdt_mlp_stack_min_intercept))
                    warn_parse("auau_nontight_bdt_mlp_stack_min_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_bdt_mlp_stack_min_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_bdt_mlp_stack_min_slope))
                    warn_parse("auau_nontight_bdt_mlp_stack_min_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_bdt_mlp_stack_max_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_bdt_mlp_stack_max_intercept))
                    warn_parse("auau_nontight_bdt_mlp_stack_max_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_bdt_mlp_stack_max_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_bdt_mlp_stack_max_slope))
                    warn_parse("auau_nontight_bdt_mlp_stack_max_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_bdt_mlp_stack_working_point_entries"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_bdt_mlp_stack_working_point_entries);
            }
            else if (StartsWithKey(line, "auau_tight_logreg_model_file"))
            {
                cfg.auau_tight_logreg_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_tight_logreg_min_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_logreg_min_intercept))
                    warn_parse("auau_tight_logreg_min_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_logreg_min_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_logreg_min_slope))
                    warn_parse("auau_tight_logreg_min_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_logreg_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_logreg_max))
                    warn_parse("auau_tight_logreg_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_logreg_apply_pt_min"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_logreg_apply_pt_min))
                    warn_parse("auau_tight_logreg_apply_pt_min", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_logreg_apply_pt_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_tight_logreg_apply_pt_max))
                    warn_parse("auau_tight_logreg_apply_pt_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_logreg_min_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_logreg_min_intercept))
                    warn_parse("auau_nontight_logreg_min_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_logreg_min_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_logreg_min_slope))
                    warn_parse("auau_nontight_logreg_min_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_logreg_max_intercept"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_logreg_max_intercept))
                    warn_parse("auau_nontight_logreg_max_intercept", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_nontight_logreg_max_slope"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_nontight_logreg_max_slope))
                    warn_parse("auau_nontight_logreg_max_slope", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_tight_logreg_working_point_entries"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.auau_tight_logreg_working_point_entries);
            }
            else if (StartsWithKey(line, "auau_bdt_training_tree"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.auau_bdt_training_tree))
                    warn_parse("auau_bdt_training_tree", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "auau_bdt_training_tree_max_entries"))
            {
                const std::string rhs = AfterColon(line);
                double val = 0.0;
                if (ParseDouble(rhs, val)) cfg.auau_bdt_training_tree_max_entries = static_cast<long long>(val);
                else warn_parse("auau_bdt_training_tree_max_entries", rhs, "expected an integer");
            }
            else if (StartsWithKey(line, "auau_bdt_npb_data_tagging"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.auau_bdt_npb_data_tagging))
                    warn_parse("auau_bdt_npb_data_tagging", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "pp_photonid_extract_only"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.pp_photonid_extract_only))
                    warn_parse("pp_photonid_extract_only", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "pp_photonid_training_tree"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.pp_photonid_training_tree))
                    warn_parse("pp_photonid_training_tree", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "pp_photonid_ppg12_filter"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.pp_photonid_ppg12_filter))
                    warn_parse("pp_photonid_ppg12_filter", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "pp_photonid_require_preselection"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.pp_photonid_require_preselection))
                    warn_parse("pp_photonid_require_preselection", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "pp_photonid_training_tree_max_entries"))
            {
                const std::string rhs = AfterColon(line);
                double val = 0.0;
                if (ParseDouble(rhs, val)) cfg.pp_photonid_training_tree_max_entries = static_cast<long long>(val);
                else warn_parse("pp_photonid_training_tree_max_entries", rhs, "expected an integer");
            }
            else if (StartsWithKey(line, "pp_photonid_source_role"))
            {
                cfg.pp_photonid_source_role = StripQuotes(AfterColon(line));
            }
            else if (StartsWithKey(line, "auau_npb_tag_delta_t_cut"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_npb_tag_delta_t_cut))
                    warn_parse("auau_npb_tag_delta_t_cut", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_npb_tag_weta_min"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_npb_tag_weta_min))
                    warn_parse("auau_npb_tag_weta_min", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_npb_tag_away_jet_pt_min"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_npb_tag_away_jet_pt_min))
                    warn_parse("auau_npb_tag_away_jet_pt_min", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_npb_tag_away_jet_dphi_min"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_npb_tag_away_jet_dphi_min))
                    warn_parse("auau_npb_tag_away_jet_dphi_min", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_npb_tag_time_sample_ns"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_npb_tag_time_sample_ns))
                    warn_parse("auau_npb_tag_time_sample_ns", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "auau_npb_mbd_t0_offset"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.auau_npb_mbd_t0_offset))
                    warn_parse("auau_npb_mbd_t0_offset", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "jet_ml_training_tree"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.jet_ml_training_tree))
                    warn_parse("jet_ml_training_tree", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "jet_ml_training_tree_max_entries"))
            {
                const std::string rhs = AfterColon(line);
                double val = 0.0;
                if (ParseDouble(rhs, val)) cfg.jet_ml_training_tree_max_entries = static_cast<long long>(val);
                else warn_parse("jet_ml_training_tree_max_entries", rhs, "expected an integer");
            }
            else if (StartsWithKey(line, "jet_ml_correction_enabled"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.jet_ml_correction_enabled))
                    warn_parse("jet_ml_correction_enabled", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "jet_ml_use_as_nominal"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.jet_ml_use_as_nominal))
                    warn_parse("jet_ml_use_as_nominal", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "jet_ml_model_file"))
            {
                cfg.jet_ml_model_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "jet_ml_features"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListStrings(rhs, cfg.jet_ml_features);
                if (cfg.jet_ml_features.empty())
                    warn_parse("jet_ml_features", rhs, "expected an inline list of feature names");
            }
            else if (StartsWithKey(line, "doPi0Analysis"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.doPi0Analysis))
                    warn_parse("doPi0Analysis", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "event_display_tree"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.event_display_tree))
                    warn_parse("event_display_tree", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "event_display_tree_max_per_bin"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseInt(rhs, cfg.event_display_tree_max_per_bin))
                {
                    warn_parse("event_display_tree_max_per_bin", rhs, "expected an integer");
                }
                if (cfg.event_display_tree_max_per_bin < 0) cfg.event_display_tree_max_per_bin = 0;
            }
            else if (StartsWithKey(line, "matching"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("pho_dr_max")) cfg.pho_dr_max = m["pho_dr_max"];
                if (m.count("jet_dr_max")) cfg.jet_dr_max = m["jet_dr_max"];
            }
            else if (StartsWithKey(line, "isolation_wp"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("aGeV"))       cfg.isoA        = m["aGeV"];
                if (m.count("bPerGeV"))    cfg.isoB        = m["bPerGeV"];
                if (m.count("sideGapGeV")) cfg.isoGap      = m["sideGapGeV"];
                if (m.count("fixedGeV"))   cfg.isoFixed    = m["fixedGeV"];
                if (m.count("truthIsoGeV"))cfg.truthIsoGeV = m["truthIsoGeV"];
                if (m.count("coneR"))      cfg.isoConeR    = m["coneR"];
                if (m.count("towerMin"))   cfg.isoTowMin   = m["towerMin"];
            }
            else if (StartsWithKey(line, "isSlidingIso"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.isSlidingIso))
                    warn_parse("isSlidingIso", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "isSlidingAndFixed"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.isSlidingAndFixed))
                    warn_parse("isSlidingAndFixed", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "fixedGeV"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.isoFixed))
                {
                    std::vector<double> v;
                    ParseInlineListDoubles(rhs, v);
                    if (!v.empty())
                    {
                        cfg.isoFixed = v.front();
                        std::ostringstream oss;
                        oss << "[CFG] fixedGeV is a list (n=" << v.size() << "); using first value = " << cfg.isoFixed;
                        info_parse(oss.str());
                    }
                    else
                    {
                        warn_parse("fixedGeV", rhs, "expected a scalar double or an inline list [..]");
                    }
                }
            }
            else if (StartsWithKey(line, "pp_iso_wp_r30"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("aGeV") && m.count("bPerGeV"))
                {
                    cfg.ppIsoWPR30.aGeV = m["aGeV"];
                    cfg.ppIsoWPR30.bPerGeV = m["bPerGeV"];
                    cfg.ppIsoWPR30.sideGapGeV = (m.count("sideGapGeV") ? m["sideGapGeV"] : cfg.isoGap);
                    cfg.ppIsoWPR30.configured = true;
                }
                else
                {
                    warn_parse("pp_iso_wp_r30", AfterColon(line), "expected inline map with aGeV and bPerGeV");
                }
            }
            else if (StartsWithKey(line, "pp_iso_wp_r40"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("aGeV") && m.count("bPerGeV"))
                {
                    cfg.ppIsoWPR40.aGeV = m["aGeV"];
                    cfg.ppIsoWPR40.bPerGeV = m["bPerGeV"];
                    cfg.ppIsoWPR40.sideGapGeV = (m.count("sideGapGeV") ? m["sideGapGeV"] : cfg.isoGap);
                    cfg.ppIsoWPR40.configured = true;
                }
                else
                {
                    warn_parse("pp_iso_wp_r40", AfterColon(line), "expected inline map with aGeV and bPerGeV");
                }
            }
            else if (StartsWithKey(line, "auau_cent_iso_wp_r30"))
            {
                parse_cent_iso_wp_list("auau_cent_iso_wp_r30", cfg.auauCentIsoWPR30);
            }
            else if (StartsWithKey(line, "auau_cent_iso_wp_r40"))
            {
                parse_cent_iso_wp_list("auau_cent_iso_wp_r40", cfg.auauCentIsoWPR40);
            }
            else if (StartsWithKey(line, "auau_cent_iso_wp"))
            {
                parse_cent_iso_wp_list("auau_cent_iso_wp", cfg.auauCentIsoWP);
            }
            else if (StartsWithKey(line, "photon_id_pre"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("e11e33_max")) cfg.pre_e11e33_max = m["e11e33_max"];
                if (m.count("et1_min"))    cfg.pre_et1_min    = m["et1_min"];
                if (m.count("et1_max"))    cfg.pre_et1_max    = m["et1_max"];
                if (m.count("e32e35_min")) cfg.pre_e32e35_min = m["e32e35_min"];
                if (m.count("e32e35_max")) cfg.pre_e32e35_max = m["e32e35_max"];
                if (m.count("weta_max"))   cfg.pre_weta_max   = m["weta_max"];
            }
            else if (StartsWithKey(line, "photon_id_tight"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("w_lo"))            cfg.tight_w_lo           = m["w_lo"];
                if (m.count("w_hi_intercept"))  cfg.tight_w_hi_intercept = m["w_hi_intercept"];
                if (m.count("w_hi_slope"))      cfg.tight_w_hi_slope     = m["w_hi_slope"];
                
                if (m.count("e11e33_min")) cfg.tight_e11e33_min = m["e11e33_min"];
                if (m.count("e11e33_max")) cfg.tight_e11e33_max = m["e11e33_max"];
                
                if (m.count("et1_min"))    cfg.tight_et1_min    = m["et1_min"];
                if (m.count("et1_max"))    cfg.tight_et1_max    = m["et1_max"];
                
                if (m.count("e32e35_min")) cfg.tight_e32e35_min = m["e32e35_min"];
                if (m.count("e32e35_max")) cfg.tight_e32e35_max = m["e32e35_max"];
            }
            else if (StartsWithKey(line, "jes3_photon_pt_bins"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListDoubles(rhs, cfg.jes3_photon_pt_bins);
                if (cfg.jes3_photon_pt_bins.size() < 2)
                    warn_parse("jes3_photon_pt_bins", rhs, "expected an inline list with >=2 edges");
            }
            else if (StartsWithKey(line, "unfold_reco_photon_pt_bins"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListDoubles(rhs, cfg.unfold_reco_photon_pt_bins);
                if (cfg.unfold_reco_photon_pt_bins.size() < 2)
                    warn_parse("unfold_reco_photon_pt_bins", rhs, "expected an inline list with >=2 edges");
            }
            else if (StartsWithKey(line, "unfold_truth_photon_pt_bins"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListDoubles(rhs, cfg.unfold_truth_photon_pt_bins);
                if (cfg.unfold_truth_photon_pt_bins.size() < 2)
                    warn_parse("unfold_truth_photon_pt_bins", rhs, "expected an inline list with >=2 edges");
            }
            else if (StartsWithKey(line, "unfold_xj_bins"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListDoubles(rhs, cfg.unfold_xj_bins);
                if (cfg.unfold_xj_bins.size() < 2)
                    warn_parse("unfold_xj_bins", rhs, "expected an inline list with >=2 edges");
            }
            else if (StartsWithKey(line, "unfold_jet_pt_binning"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("start")) cfg.unfold_jet_pt_start = m["start"];
                if (m.count("stop"))  cfg.unfold_jet_pt_stop  = m["stop"];
                if (m.count("step"))  cfg.unfold_jet_pt_step  = m["step"];
            }
            else if (StartsWithKey(line, "leading_response_family"))
            {
                cfg.leading_response_family = detail::trim(AfterColon(line));
                if (cfg.leading_response_family == "nominal")
                {
                    cfg.leading_response_family.clear();
                }
                else if (!cfg.leading_response_family.empty() &&
                         cfg.leading_response_family != "sam_compat")
                {
                    throw std::runtime_error(
                        "leading_response_family must be empty, nominal, or sam_compat");
                }
            }
            else if (StartsWithKey(line, "require_towerinfo_truth_matching"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.require_towerinfo_truth_matching))
                    warn_parse("require_towerinfo_truth_matching", rhs, "expected true/false");
            }
        }

        if (cfg.auauCentIsoWP.empty() && yaml_has_key("auau_cent_iso_wp"))
            recover_cent_iso_wp_block("auau_cent_iso_wp", cfg.auauCentIsoWP);
        if (cfg.auauCentIsoWPR30.empty() && yaml_has_key("auau_cent_iso_wp_r30"))
            recover_cent_iso_wp_block("auau_cent_iso_wp_r30", cfg.auauCentIsoWPR30);
        if (cfg.auauCentIsoWPR40.empty() && yaml_has_key("auau_cent_iso_wp_r40"))
            recover_cent_iso_wp_block("auau_cent_iso_wp_r40", cfg.auauCentIsoWPR40);

        auto require_loaded_cent_iso_wp = [&](const std::string& key, const std::vector<Config::CentIsoWP>& values)
        {
            if (!yaml_has_key(key) || !values.empty()) return;
            std::ostringstream oss;
            oss << "YAML contains '" << key
                << "' but no centrality-dependent isolation entries were loaded; "
                << "refusing to run with fixed-isolation fallback.";
            throw std::runtime_error(oss.str());
        };
        require_loaded_cent_iso_wp("auau_cent_iso_wp", cfg.auauCentIsoWP);
        require_loaded_cent_iso_wp("auau_cent_iso_wp_r30", cfg.auauCentIsoWPR30);
        require_loaded_cent_iso_wp("auau_cent_iso_wp_r40", cfg.auauCentIsoWPR40);
        
        if (yamlV > 0)
        {
            std::cout << "[CFG] YAML parsing completed OK: " << cfg.yamlPath
            << (strict ? " (strict)" : "") << "\n";
        }
        
        return cfg;
    }
    
    inline void ExpandUniformEdges(std::vector<double>& edges, double start, double stop, double step)
    {
        edges.clear();
        if (!(std::isfinite(start) && std::isfinite(stop) && std::isfinite(step))) return;
        if (step <= 0.0) return;
        if (stop <= start) return;
        
        const int n = (int) std::llround((stop - start) / step);
        edges.reserve((std::size_t)n + 2);
        for (int i = 0; i <= n; ++i)
        {
            edges.push_back(start + step * (double)i);
        }
        if (edges.empty() || std::fabs(edges.back() - stop) > 1e-9) edges.push_back(stop);
    }
}

namespace idfanout
{
    struct Entry
    {
        std::string outRoot;
        std::string cfgTag;
        std::string preselection;
        std::string tight;
        std::string nonTight;
        bool hasIsoOverride = false;
        double coneR = 0.0;
        bool isSlidingIso = false;
        double fixedGeV = 0.0;
    };

    inline std::vector<std::string> SplitPipe(const std::string& line)
    {
        std::vector<std::string> out;
        std::string cur;
        std::istringstream ss(line);
        while (std::getline(ss, cur, '|')) out.push_back(detail::trim(cur));
        return out;
    }

    inline std::vector<Entry> LoadFromEnv()
    {
        std::vector<Entry> entries;
        const char* rawPath = std::getenv("RJ_ID_FANOUT_FILE");
        if (!rawPath || !*rawPath) return entries;

        const std::string path = detail::trim(std::string(rawPath));
        if (path.empty()) return entries;

        std::ifstream in(path);
        if (!in.is_open())
        {
            detail::bail("RJ_ID_FANOUT_FILE is set but cannot be opened: " + path);
        }

        for (std::string line; std::getline(in, line); )
        {
            line = detail::trim(line);
            if (line.empty() || line[0] == '#') continue;
            const auto cols = SplitPipe(line);
            if (!(cols.size() == 5 || cols.size() == 8))
            {
                detail::bail("Bad RJ_ID_FANOUT_FILE row. Expected outRoot|cfgTag|preselection|tight|nonTight or outRoot|cfgTag|preselection|tight|nonTight|coneR|isSlidingIso|fixedGeV, got: " + line);
            }

            Entry e;
            e.outRoot = cols[0];
            e.cfgTag = cols[1];
            e.preselection = yamlcfg::NormalizePreselectionMode(cols[2]);
            e.tight = yamlcfg::NormalizeTightMode(cols[3]);
            e.nonTight = yamlcfg::NormalizeNonTightMode(cols[4]);

            if (e.outRoot.empty()) detail::bail("RJ_ID_FANOUT_FILE row has empty outRoot: " + line);
            if (!yamlcfg::IsPreselectionMode(e.preselection))
                detail::bail("RJ_ID_FANOUT_FILE row has invalid preselection: " + cols[2]);
            if (!yamlcfg::IsTightMode(e.tight))
                detail::bail("RJ_ID_FANOUT_FILE row has invalid tight: " + cols[3]);
            if (!yamlcfg::IsNonTightMode(e.nonTight))
                detail::bail("RJ_ID_FANOUT_FILE row has invalid nonTight: " + cols[4]);

            if (cols.size() == 8)
            {
                try
                {
                    e.coneR = std::stod(cols[5]);
                    e.fixedGeV = std::stod(cols[7]);
                }
                catch (...)
                {
                    detail::bail("RJ_ID_FANOUT_FILE row has invalid coneR/fixedGeV numeric fields: " + line);
                }
                if (!yamlcfg::ParseBool(cols[6], e.isSlidingIso))
                    detail::bail("RJ_ID_FANOUT_FILE row has invalid isSlidingIso field: " + cols[6]);
                if (!std::isfinite(e.coneR) || !std::isfinite(e.fixedGeV))
                    detail::bail("RJ_ID_FANOUT_FILE row has non-finite coneR/fixedGeV: " + line);
                e.hasIsoOverride = true;
            }

            entries.push_back(e);
        }

        if (entries.empty())
        {
            detail::bail("RJ_ID_FANOUT_FILE was provided but contained no usable rows: " + path);
        }
        return entries;
    }

    inline void ReplaceOrAppendScalar(std::string& text,
                                      const std::string& key,
                                      const std::string& value)
    {
        std::istringstream in(text);
        std::ostringstream out;
        bool replaced = false;
        for (std::string line; std::getline(in, line); )
        {
            std::string trimmed = detail::trim(line);
            if (!replaced && trimmed.rfind(key + ":", 0) == 0)
            {
                out << key << ": " << value << '\n';
                replaced = true;
            }
            else
            {
                out << line << '\n';
            }
        }
        if (!replaced) out << key << ": " << value << '\n';
        text = out.str();
    }

    inline std::string YAMLForEntry(const std::string& baseYaml, const Entry& e)
    {
        std::string text = baseYaml;
        ReplaceOrAppendScalar(text, "preselection", e.preselection);
        ReplaceOrAppendScalar(text, "tight", e.tight);
        ReplaceOrAppendScalar(text, "nonTight", e.nonTight);
        if (e.hasIsoOverride)
        {
            ReplaceOrAppendScalar(text, "coneR", detail::fmt(e.coneR, 3));
            ReplaceOrAppendScalar(text, "isSlidingIso", e.isSlidingIso ? "true" : "false");
            ReplaceOrAppendScalar(text, "isSlidingAndFixed", "false");
            ReplaceOrAppendScalar(text, "fixedGeV", detail::fmt(e.fixedGeV, 3));
        }
        if (!e.cfgTag.empty()) ReplaceOrAppendScalar(text, "analysis_cfg_tag", e.cfgTag);
        return text;
    }
}

// ============================================================================
// Ensure the composite nodes JetCalib hard-requires exist.
// JetCalib::CreateNodeTree requires BOTH:
//   - PHCompositeNode "ANTIKT"
//   - PHCompositeNode "TOWER"
// Your pp DSTs often don't have "TOWER", so JetCalib aborts the run.
// ============================================================================
class EnsureJetCalibNodes final : public SubsysReco
{
 public:
  explicit EnsureJetCalibNodes(const std::string& name="EnsureJetCalibNodes")
    : SubsysReco(name) {}

  int InitRun(PHCompositeNode* topNode) override
  {
    PHNodeIterator iter(topNode);

    auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      std::cerr << "[EnsureJetCalibNodes] FATAL: DST node not found\n";
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // JetCalib requires ANTIKT to exist
    auto* antiktNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "ANTIKT"));
    if (!antiktNode)
    {
      antiktNode = new PHCompositeNode("ANTIKT");
      dstNode->addNode(antiktNode);
      if (Verbosity() > 0)
        std::cout << "[EnsureJetCalibNodes] created DST/ANTIKT composite node\n";
    }

    // JetCalib requires TOWER to exist (this is what you're missing)
    auto* towerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "TOWER"));
    if (!towerNode)
    {
      towerNode = new PHCompositeNode("TOWER");
      dstNode->addNode(towerNode);
      if (Verbosity() > 0)
        std::cout << "[EnsureJetCalibNodes] created DST/TOWER composite node\n";
    }

    return Fun4AllReturnCodes::EVENT_OK;
  }
};

class ProcessEnvSetter final : public SubsysReco
{
 public:
  ProcessEnvSetter(const std::string& name,
                   const std::string& key,
                   const std::string& value)
    : SubsysReco(name)
    , m_key(key)
    , m_value(value)
  {}

  int Init(PHCompositeNode* /*topNode*/) override
  {
    setenv(m_key.c_str(), m_value.c_str(), 1);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int InitRun(PHCompositeNode* /*topNode*/) override
  {
    setenv(m_key.c_str(), m_value.c_str(), 1);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int process_event(PHCompositeNode* /*topNode*/) override
  {
    setenv(m_key.c_str(), m_value.c_str(), 1);
    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  std::string m_key;
  std::string m_value;
};


class TowerAudit final : public SubsysReco
{
 public:
  explicit TowerAudit(const std::string& name,
                      const std::string& nodeCEMC,
                      const std::string& nodeHCIN,
                      const std::string& nodeHCOUT,
                      int maxPrintEvents = 5)
    : SubsysReco(name)
    , m_nodeCEMC(nodeCEMC)
    , m_nodeHCIN(nodeHCIN)
    , m_nodeHCOUT(nodeHCOUT)
    , m_maxPrint(maxPrintEvents)
  {}

  int InitRun(PHCompositeNode* topNode) override
  {
    // Just verify nodes exist at InitRun (they may still be filled per-event)
    auto* cemc = findNode::getClass<TowerInfoContainer>(topNode, m_nodeCEMC);
    auto* hcin = findNode::getClass<TowerInfoContainer>(topNode, m_nodeHCIN);
    auto* hcot = findNode::getClass<TowerInfoContainer>(topNode, m_nodeHCOUT);

    if (Verbosity() > 0)
    {
      std::cout << "[" << Name() << "::InitRun]"
                << " nodes:"
                << " CEMC=" << m_nodeCEMC << (cemc ? "(OK)" : "(MISSING)")
                << " HCALIN=" << m_nodeHCIN << (hcin ? "(OK)" : "(MISSING)")
                << " HCALOUT=" << m_nodeHCOUT << (hcot ? "(OK)" : "(MISSING)")
                << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int process_event(PHCompositeNode* topNode) override
  {
    ++m_evt;
    if (m_evt > m_maxPrint) return Fun4AllReturnCodes::EVENT_OK;

    auto auditOne = [&](const char* label, const std::string& node)
    {
      auto* cont = findNode::getClass<TowerInfoContainer>(topNode, node);
      if (!cont)
      {
        std::cout << "[" << Name() << "] evt=" << m_evt
                  << " " << label << " node=" << node << " MISSING\n";
        return;
      }

      const unsigned int nt = cont->size();
      unsigned int nNull = 0, nGood = 0, nBad = 0;
      unsigned int nEpos = 0, nEposGood = 0, nEposBad = 0;

      double sumE = 0.0, sumEgood = 0.0, sumEbad = 0.0;
      double sumT = 0.0, sumTgood = 0.0, sumTbad = 0.0;
      unsigned int nT = 0, nTgood = 0, nTbad = 0;

      // sample the whole container (24576 for CEMC is fine for a few events)
      for (unsigned int ch = 0; ch < nt; ++ch)
      {
        TowerInfo* t = cont->get_tower_at_channel(ch);
        if (!t) { ++nNull; continue; }

        const bool good = t->get_isGood();
        const float e   = t->get_energy();
        const float tim = t->get_time();

        if (good) ++nGood; else ++nBad;

        if (std::isfinite(e))
        {
          sumE += e;
          if (e > 0) ++nEpos;
          if (good) { sumEgood += e; if (e > 0) ++nEposGood; }
          else      { sumEbad  += e; if (e > 0) ++nEposBad;  }
        }

        if (std::isfinite(tim))
        {
          ++nT;
          sumT += tim;
          if (good) { ++nTgood; sumTgood += tim; }
          else      { ++nTbad;  sumTbad  += tim; }
        }
      }

      const double fracBad = (nt > 0) ? (double)nBad / (double)nt : 0.0;
      const double meanT   = (nT > 0) ? sumT / (double)nT : 0.0;
      const double meanTg  = (nTgood > 0) ? sumTgood / (double)nTgood : 0.0;
      const double meanTb  = (nTbad > 0) ? sumTbad / (double)nTbad : 0.0;

      std::cout << "[" << Name() << "] evt=" << m_evt << " " << label
                << " node=" << node
                << " ntowers=" << nt
                << " null=" << nNull
                << " good=" << nGood
                << " bad=" << nBad
                << " fracBad=" << std::fixed << std::setprecision(4) << fracBad
                << " | sumE=" << std::setprecision(3) << sumE
                << " (good=" << sumEgood << ", bad=" << sumEbad << ")"
                << " | E>0: all=" << nEpos << " good=" << nEposGood << " bad=" << nEposBad
                << " | meanTime=" << std::setprecision(3) << meanT
                << " (good=" << meanTg << ", bad=" << meanTb << ")"
                << "\n";
    };

    auditOne("CEMC",   m_nodeCEMC);
    auditOne("HCALIN", m_nodeHCIN);
    auditOne("HCALOUT",m_nodeHCOUT);

    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  std::string m_nodeCEMC;
  std::string m_nodeHCIN;
  std::string m_nodeHCOUT;
  int m_maxPrint = 5;
  int m_evt = 0;
};



class JetCalibOneEventProbe final : public SubsysReco
{
 public:
  JetCalibOneEventProbe(const std::string& name,
                        const std::string& rawNode,
                        const std::string& calibNode,
                        int maxJetsToPrint = 12)
    : SubsysReco(name)
    , m_rawNode(rawNode)
    , m_calibNode(calibNode)
    , m_maxJets(maxJetsToPrint)
  {}

  int process_event(PHCompositeNode* topNode) override
  {
    // Hard gate: ONLY show anything at ultra-verbose
    if (Verbosity() < 20) return Fun4AllReturnCodes::EVENT_OK;

    // Print only once for the whole job
    if (m_fired) return Fun4AllReturnCodes::EVENT_OK;

    auto* raw   = findNode::getClass<JetContainer>(topNode, m_rawNode);
    auto* calib = findNode::getClass<JetContainer>(topNode, m_calibNode);

      // Only trigger when both exist and actually have jets
      if (!raw || !calib) return Fun4AllReturnCodes::EVENT_OK;
      if (raw->size() == 0 || calib->size() == 0) return Fun4AllReturnCodes::EVENT_OK;

      // Grab z-vertex (MBD) JetCalib uses; REQUIRE |zvtx| < 30 cm
      float zvtx = -999.0f;
      auto* vtxmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
      if (vtxmap && !vtxmap->empty())
      {
        GlobalVertex* vtx = vtxmap->begin()->second;
        if (vtx)
        {
          auto mbdStart = vtx->find_vertexes(GlobalVertex::MBD);
          auto mbdEnd   = vtx->end_vertexes();
          for (auto it = mbdStart; it != mbdEnd; ++it)
          {
            const auto& [type, vec] = *it;
            if (type != GlobalVertex::MBD) continue;
            for (const auto* vv : vec)
            {
              if (!vv) continue;
              zvtx = vv->get_z();
            }
          }
        }
      }

      // Require a well-defined vertex AND a central-ish event
      if (!std::isfinite(zvtx)) return Fun4AllReturnCodes::EVENT_OK;
      if (std::fabs(zvtx) >= 30.0f) return Fun4AllReturnCodes::EVENT_OK;

      const int nRaw   = (int) raw->size();
      const int nCalib = (int) calib->size();
      const int nScan  = std::min(nRaw, nCalib);

      // Only print jets with pt_raw >= 5 GeV
      const float ptMinPrint = 5.0f;

      // First pass: count how many jets pass pt threshold (to decide whether to print at all)
      int nPass = 0;
      for (int i = 0; i < nScan; ++i)
      {
        const Jet* jr = raw->get_jet(i);
        if (!jr) continue;
        if (jr->get_pt() >= ptMinPrint) ++nPass;
      }

      // If nothing interesting, don't fire and keep searching future events
      if (nPass == 0) return Fun4AllReturnCodes::EVENT_OK;

      // ---------------- ANSI helpers ----------------
      const char* RST  = "\033[0m";
      const char* BOLD = "\033[1m";
      const char* DIM  = "\033[2m";

      const char* C_RAW   = "\033[38;5;45m";   // cyan-ish
      const char* C_CAL   = "\033[38;5;208m";  // orange-ish
      const char* C_INFO  = "\033[38;5;111m";  // light blue
      const char* C_WARN  = "\033[38;5;220m";  // yellow
      const char* C_GOOD  = "\033[38;5;82m";   // green
      const char* C_BAD   = "\033[38;5;196m";  // red
      const char* C_BAR   = "\033[38;5;244m";  // gray

      auto scaleColor = [&](float sc) -> const char*
      {
        if (!std::isfinite(sc) || sc < 0.0f) return C_WARN;
        if (sc > 1.05f) return C_GOOD;
        if (sc < 0.95f) return C_BAD;
        return C_WARN;
      };

      std::cout << "\n" << C_BAR << "====================================================================" << RST << "\n";
      std::cout << BOLD << C_INFO << "[JES PROBE] One-event before/after JetCalib (filtered)" << RST << "\n";

      std::cout << "  " << BOLD << C_RAW << "RAW" << RST
                << "   node: " << C_RAW << m_rawNode << RST << "  " << DIM << "(n=" << nRaw << ")" << RST << "\n";
      std::cout << "  " << BOLD << C_CAL << "CAL" << RST
                << "   node: " << C_CAL << m_calibNode << RST << "  " << DIM << "(n=" << nCalib << ")" << RST << "\n";

      std::cout << "  " << BOLD << "zvtx(MBD):" << RST << " " << C_INFO
                << std::fixed << std::setprecision(2) << zvtx << " cm" << RST
                << "  " << DIM << "(|z|<30 required)" << RST << "\n";

      std::cout << "  " << BOLD << "pt_raw >= " << RST << C_WARN
                << std::fixed << std::setprecision(1) << ptMinPrint << " GeV" << RST
                << "  " << DIM << "(printing " << nPass << " jets; maxRows=" << m_maxJets << ")" << RST << "\n";

      std::cout << C_BAR << "--------------------------------------------------------------------" << RST << "\n";

      // -------- fixed-width table layout (headers match rows) ----------
      const int W_I   = 3;
      const int W_PT  = 8;
      const int W_ETA = 8;
      const int W_PHI = 8;
      const int W_SC  = 6;

      // Header row (use SAME widths as data)
      std::cout
        << "  " << BOLD << std::setw(W_I) << "i" << RST << " " << C_BAR << "|" << RST << " "
        << BOLD << C_RAW << std::setw(W_PT)  << "pt_raw"  << RST << " "
        << BOLD << C_RAW << std::setw(W_ETA) << "eta_raw" << RST << " "
        << BOLD << C_RAW << std::setw(W_PHI) << "phi_raw" << RST << "  "
        << C_BAR << "=>" << RST << "  "
        << BOLD << C_CAL << std::setw(W_PT)  << "pt_cal"  << RST << " "
        << BOLD << C_CAL << std::setw(W_ETA) << "eta_cal" << RST << " "
        << BOLD << C_CAL << std::setw(W_PHI) << "phi_cal" << RST << "  "
        << C_BAR << "||" << RST << " "
        << BOLD << std::setw(W_SC) << "scale" << RST
        << "\n";

      // Separator line that matches the table width exactly
      auto rep = [](int n, char c) { return std::string(std::max(0, n), c); };
      const int tableWidth =
          2  // leading two spaces
        + W_I + 1 + 1 + 1  // i + space + '|' + space
        + W_PT + 1 + W_ETA + 1 + W_PHI  // raw cols with spaces
        + 2 + 2 + 2  // two spaces + "=>" + two spaces (approx; we print literal)
        + W_PT + 1 + W_ETA + 1 + W_PHI  // cal cols
        + 2 + 2 + 1 + W_SC; // two spaces + "||" + space + scale

      std::cout << C_BAR << rep(tableWidth, '-') << RST << "\n";

      int nPrinted = 0;
      for (int i = 0; i < nScan; ++i)
      {
        const Jet* jr = raw->get_jet(i);
        const Jet* jc = calib->get_jet(i);
        if (!jr || !jc) continue;

        const float pr = jr->get_pt();
        if (pr < ptMinPrint) continue;

        const float pc = jc->get_pt();
        const float sc = (pr > 0.0f) ? (pc / pr) : -1.0f;

        std::cout
          << "  " << BOLD << std::setw(W_I) << i << RST << " " << C_BAR << "|" << RST << " "
          << C_RAW
          << std::setw(W_PT)  << std::fixed << std::setprecision(3) << pr << " "
          << std::setw(W_ETA) << std::setprecision(3) << jr->get_eta() << " "
          << std::setw(W_PHI) << std::setprecision(3) << jr->get_phi()
          << RST << "  "
          << C_BAR << "=>" << RST << "  "
          << C_CAL
          << std::setw(W_PT)  << std::setprecision(3) << pc << " "
          << std::setw(W_ETA) << std::setprecision(3) << jc->get_eta() << " "
          << std::setw(W_PHI) << std::setprecision(3) << jc->get_phi()
          << RST << "  "
          << C_BAR << "||" << RST << " "
          << scaleColor(sc) << BOLD << std::setw(W_SC) << std::setprecision(3) << sc << RST
          << "\n";

        ++nPrinted;
        if (nPrinted >= m_maxJets) break;
      }

      std::cout << C_BAR << "====================================================================" << RST << "\n\n";

      m_fired = true;
      return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  std::string m_rawNode;
  std::string m_calibNode;
  int m_maxJets = 12;
  bool m_fired = false;
};



//======================================================================
//  The actual steering macro
//======================================================================
void Fun4All_recoilJets_unified_impl(const int   nEvents   =  0,
                                     const char* listFile  = "input_files.list",
                                     const char* outRoot   = "TrigPlot.root",
                                     const bool  verbose   = false)
{
    //--------------------------------------------------------------------
    // 0.  Banner & basic environment sanity
    //--------------------------------------------------------------------
    if (verbose) {
        std::cout << "\n>>> Fun4All_recoilJets – ana.495 driver <<<\n"
        << "    Input list : " << listFile  << '\n'
        << "    Output file: " << outRoot   << '\n'
        << "    nEvents    : " << nEvents   << (nEvents==0? " (all)\n":"\n");
    }
    
    Fun4AllServer* se = Fun4AllServer::instance();
    Fun4AllInputManager* permittedRepeatingPedestalInputManager = nullptr;
    if (!se) detail::bail("unable to obtain Fun4AllServer instance!");
    
    
    auto env_lower = [](const char* key, const std::string& def = std::string{}) -> std::string
    {
        const char* raw = std::getenv(key);
        std::string s = raw ? detail::trim(std::string(raw)) : def;
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
        return s;
    };

    auto env_truthy_local = [&](const char* key) -> bool
    {
        const std::string v = env_lower(key);
        return v == "1" || v == "true" || v == "yes" || v == "on";
    };

    auto env_bool_local = [&](const char* key, bool def) -> bool
    {
        const std::string v = env_lower(key, def ? "1" : "0");
        if (v == "1" || v == "true" || v == "yes" || v == "on") return true;
        if (v == "0" || v == "false" || v == "no" || v == "off") return false;
        return def;
    };

    const char* analysisJetNodeSuffixRaw =
        std::getenv("RJ_ANALYSIS_JET_NODE_SUFFIX");
    const std::string analysisJetNodeSuffix = analysisJetNodeSuffixRaw
        ? detail::trim(std::string(analysisJetNodeSuffixRaw))
        : std::string{};
    if (!RJReplayFoundationV1::validAnalysisJetNodeSuffix(
            analysisJetNodeSuffix))
    {
        detail::bail(
            "RJ_ANALYSIS_JET_NODE_SUFFIX must be empty or an underscore-prefixed "
            "alphanumeric token no longer than 32 characters");
    }
    if (env_truthy_local("RJ_REPLAY_FOUNDATION_V1") &&
        analysisJetNodeSuffix.empty())
    {
        detail::bail(
            "RJ_REPLAY_FOUNDATION_V1 requires a nonempty "
            "RJ_ANALYSIS_JET_NODE_SUFFIX so reconstructed jets cannot append "
            "into input-DST nodes");
    }

    // A bounded shower-contract diagnostic needs the ordinary Au+Au
    // calorimeter and photon-candidate reconstruction, but it does not need
    // final recoil jets or rebuilt truth jets.  Keep this strictly
    // environment-gated so production steering remains unchanged.
    const bool auauCandidateSkimOnly =
        env_bool_local("RJ_AUAU_CANDIDATE_SKIM_ONLY", false);

    //--------------------------------------------------------------------
    // 1.  Parse the file list & determine run / segment
    //--------------------------------------------------------------------
    std::ifstream list(listFile);
    if (!list.is_open())
        detail::bail("cannot open input list \"" + std::string(listFile) + "\"");
    
    // ---------------------------------------------------------------
    // Parse list file
    //   DATA:  1 column  -> calo DST
    //          2 columns -> pp PPG12: DST_Jet + DST_JETCALO
    //                    -> AuAu data policy: DST_JET + DST_JETCALO
    //                    -> AuAu MB gate: calo DST + DST_ZDC_RAW
    //   SIM :  5 columns -> calo + G4Hits + (truth jets) + global + mbd_epd
    //
    // NOTE:
    //   The second token is interpreted by mode:
    //     pp DATA: DST_JETCALO, with column 1 carrying DST_Jet
    //     AuAu DATA: DST_JETCALO for the paired jet stream, or DST_ZDC_RAW
    //                 for the legacy MB-gate stream
    //     SIM : G4Hits
    // ---------------------------------------------------------------
    std::vector<std::string> filesCalo;
    std::vector<std::string> filesZdc;
    std::vector<std::string> filesG4;
    std::vector<std::string> filesJets;
    std::vector<std::string> filesGlobal;
    std::vector<std::string> filesMbd;
    
    for (std::string line; std::getline(list, line); )
    {
        line = detail::trim(line);
        if (line.empty()) continue;
        if (!line.empty() && line[0] == '#') continue;
        
        // Columns:
        //   DATA : <DST_CALO_CLUSTER> [<DST_ZDC_RAW>]
        //      or <DST_Jet> <DST_JETCALO> for PPG12 pp-data parity
        //      or <DST_JET> <DST_JETCALO> for the AuAu paired data policy
        //   SIM  : <DST_CALO_CLUSTER> <G4Hits> <DST_JETS> <DST_GLOBAL> <DST_MBD_EPD>
        std::istringstream iss(line);
        std::string fCalo, fAux1, fJets, fGlobal, fMbd;
        iss >> fCalo >> fAux1 >> fJets >> fGlobal >> fMbd;
        
        if (fCalo.empty()) continue;
        const std::string fCaloNorm = (fCalo != "NONE") ? fCalo : std::string{};

        filesCalo.emplace_back(fCaloNorm);
        
        // Keep vectors key-aligned (same length as filesCalo):
        // empty string == "not provided on this line"
        filesZdc.emplace_back((!fAux1.empty() && fAux1 != "NONE") ? fAux1 : std::string{});
        filesG4.emplace_back((!fAux1.empty() && fAux1 != "NONE") ? fAux1 : std::string{});
        filesJets.emplace_back((!fJets.empty() && fJets != "NONE") ? fJets : std::string{});
        filesGlobal.emplace_back((!fGlobal.empty() && fGlobal != "NONE") ? fGlobal : std::string{});
        filesMbd.emplace_back((!fMbd.empty() && fMbd != "NONE") ? fMbd : std::string{});
    }
    
    if (filesCalo.empty())
        detail::bail("input list \"" + std::string(listFile) + "\" is empty");
    
    std::string firstFile; // used for GetRunSegment
    auto captureFirstInput = [&](const std::vector<std::string>& paths)
    {
        if (!firstFile.empty()) return;
        for (const auto& path : paths)
        {
            if (!path.empty() && path != "NONE")
            {
                firstFile = path;
                return;
            }
        }
    };
    captureFirstInput(filesCalo);
    captureFirstInput(filesG4);
    captureFirstInput(filesJets);
    captureFirstInput(filesGlobal);
    captureFirstInput(filesMbd);
    captureFirstInput(filesZdc);
    if (firstFile.empty())
        detail::bail("input list \"" + std::string(listFile) + "\" contains no usable file tokens");
    
    // Dataset / input-mode detection
    //  - RJ_DATASET drives analysis mode: isPP | isPPrun25 | isAuAu | isSim | isSimEmbedded
    //  - RJ_CALO_INPUT_MODE optionally overrides input provenance: jetcalo | calofitting | simdst
    bool isSim = false;
    bool isSimEmbedded = false;
    bool isPPrun25 = false;
    bool isAuAuRequested = false;
    std::string datasetToken = env_lower("RJ_DATASET", "ispp");
    if (datasetToken == "issimembedded" || datasetToken == "simembedded" ||
        datasetToken == "issimembeddedinclusive" || datasetToken == "simembeddedinclusive")
    {
        isSim = true;
        isSimEmbedded = true;
        isAuAuRequested = true;
    }
    else if (datasetToken == "issim" || datasetToken == "sim"
             || datasetToken == "issimjet5" || datasetToken == "simjet5"
             || datasetToken == "issiminclusive" || datasetToken == "siminclusive"
             || datasetToken == "issimmb" || datasetToken == "simmb")
    {
        isSim = true;
    }
    else if (datasetToken == "ispprun25" || datasetToken == "pprun25" || datasetToken == "pp25")
    {
        isPPrun25 = true;
    }
    else if (datasetToken == "isauau" || datasetToken == "auau" || datasetToken == "aa")
    {
        isAuAuRequested = true;
    }
    else
    {
        datasetToken = "ispp";
    }
    
    if (!isSim)
    {
        if (const char* f = std::getenv("RJ_IS_SIM"))
            isSim = (std::atoi(f) != 0);
        if (isSim)
        {
            datasetToken = "issim";
            isPPrun25 = false;
            isAuAuRequested = false;
        }
    }
    
    if (verbose)
    {
        std::cout << "[FLOW] dataset token resolution:"
        << " RJ_DATASET=" << datasetToken
        << " | isSim=" << (isSim ? "true" : "false")
        << " | isSimEmbedded=" << (isSimEmbedded ? "true" : "false")
        << " | isPPrun25=" << (isPPrun25 ? "true" : "false")
        << " | isAuAuRequested=" << (isAuAuRequested ? "true" : "false")
        << std::endl;
    }
    
    auto detect_input_mode = [&](const std::string& firstFileLower) -> std::string
    {
        std::string forced = env_lower("RJ_CALO_INPUT_MODE");
        if (forced == "jetcalo" || forced == "calofitting" || forced == "simdst")
        {
            if (isSimEmbedded && forced == "simdst" && !env_truthy_local("RJ_ALLOW_SIMEMBEDDED_SIMDST"))
            {
                detail::bail(
                    "RJ_CALO_INPUT_MODE=simdst would skip Process_Calo_Calib for embedded DST_CALO inputs. "
                    "The Blair-style embedded AuAu contract builds TOWERINFO_CALIB_{CEMC,HCALIN,HCALOUT} "
                    "and CLUSTERINFO_CEMC from DST_CALO before RetowerCEMC. Set RJ_ALLOW_SIMEMBEDDED_SIMDST=1 "
                    "only for a foreground diagnostic proving those nodes already exist.");
            }
            return forced;
        }
        if (isSim && !isSimEmbedded) return "simdst";
        if (isSimEmbedded)
        {
            if (firstFileLower.find("dst_calo_") != std::string::npos) return "jetcalo";
            if (firstFileLower.find("calofitting") != std::string::npos) return "calofitting";
            return "jetcalo";
        }
        if (firstFileLower.find("calofitting") != std::string::npos) return "calofitting";
        if (firstFileLower.find("jetcalo") != std::string::npos) return "jetcalo";
        if (firstFileLower.find("dst_jet") != std::string::npos) return "jetcalo";
        if (isPPrun25 || isAuAuRequested) return "calofitting";
        return "jetcalo";
    };
    
    std::string firstFileLower = firstFile;
    std::transform(firstFileLower.begin(), firstFileLower.end(), firstFileLower.begin(), [](unsigned char c){ return std::tolower(c); });
    std::string caloInputMode = detect_input_mode(firstFileLower);
    
    auto runSegPair = Fun4AllUtils::GetRunSegment(firstFile);
    int run = runSegPair.first;
    int seg = runSegPair.second;
    
    // MINIMAL / SCOPED FIX:
    // Fun4AllUtils::GetRunSegment can mis-read embedded JetXX filenames like
    //   DST_CALO_Jet20-...-data-00054404-00001-00056.root
    // and return the embedded segment token as the "run".
    // For embedded/data calibration we need the real data run (>1000), so
    // repair that here whenever the parsed run is not a valid DATA timestamp.
    const std::size_t slashPos = firstFile.find_last_of("/\\");
    const std::string baseName = (slashPos == std::string::npos) ? firstFile : firstFile.substr(slashPos + 1);
    
    const bool useEmbeddedInclusiveJetRunParse =
    isSimEmbedded &&
    baseName.compare(0, std::string("DST_CALO_Jet").size(), "DST_CALO_Jet") == 0 &&
    baseName.find("-data-") != std::string::npos;
    
    if (useEmbeddedInclusiveJetRunParse && run <= 1000)
    {
        const std::string dataTag = "-data-";
        const std::size_t dataPos = baseName.find(dataTag);
        const std::string tail = baseName.substr(dataPos + dataTag.size());
        const std::size_t dash1 = tail.find('-');
        
        if (dash1 != std::string::npos)
        {
            const std::string runTok = tail.substr(0, dash1);
            const std::string rest = tail.substr(dash1 + 1);
            const std::size_t segEnd = rest.find_first_of("-.");
            const std::string segTok = rest.substr(0, segEnd);
            
            try
            {
                const int parsedRun = std::stoi(runTok);
                const int parsedSeg = std::stoi(segTok);
                if (parsedRun > 0)
                {
                    run = parsedRun;
                    seg = (parsedSeg >= 0) ? parsedSeg : 0;
                }
            }
            catch (...) { /* keep existing fallback logic below */ }
        }
    }
    
    if (run <= 0)
    {
        if (isSim)
        {
            // If the filename doesn't encode a run number, use a safe CDB timestamp.
            run = 1;
            seg = 0;
            
            // Optional override if you ever want it:
            //   export RJ_CDB_TIMESTAMP=XXXX
            if (const char* ts = std::getenv("RJ_CDB_TIMESTAMP"))
            {
                try
                {
                    const int t = std::stoi(ts);
                    if (t > 0) run = t;
                }
                catch (...) { /* keep default */ }
            }
            
            if (verbose)
                std::cout << "[INFO] Simulation mode: run/seg not in filename → using TIMESTAMP=" << run << "\n";
        }
        else
        {
            detail::bail("failed to extract run number from first file: " + firstFile);
        }
    }
    
    if (verbose)
        std::cout << "[INFO] Run=" << run << "  Seg=" << seg
        << "  (" << filesCalo.size() << " files)\n";
    
    // --------------------------------------------------------------------
    // Global verbosity control (RJ_VERBOSITY from env; defaults to 10;
    // Condor detection → 0). Also silences std::cout/cerr globally when 0.
    // --------------------------------------------------------------------
    int vlevel = 10;
    if (const char* venv = std::getenv("RJ_VERBOSITY"))
    {
        try { vlevel = std::stoi(venv); } catch (...) { vlevel = 10; }
    }
    else
    {
        // Detect Condor jobs defensively; default to quiet in batch.
        if (std::getenv("_CONDOR_SCRATCH_DIR") || std::getenv("_CONDOR_JOB_AD"))
            vlevel = 0;
    }
    
    // RAII silence for global std::cout/cerr if vlevel==0
    struct ScopedSilence {
        std::ofstream   sink;
        std::streambuf* cout_save = nullptr;
        std::streambuf* cerr_save = nullptr;
        bool active = false;
        void enable() {
            if (active) return;
            sink.open("/dev/null");
            cout_save = std::cout.rdbuf(sink.rdbuf());
            cerr_save = std::cerr.rdbuf(sink.rdbuf());
            active = true;
        }
        void disable() {
            if (!active) return;
            std::cout.rdbuf(cout_save);
            std::cerr.rdbuf(cerr_save);
            active = false;
        }
        ~ScopedSilence() {
            disable();
        }
    } _silence;
    
    if (vlevel == 0) _silence.enable();
    
    
    if (const char* fenv = std::getenv("RJ_F4A_VERBOSE"))
    {
        const int fv = std::atoi(fenv);
        if (fv > 0) se->Verbosity(fv);
    }
    
    // --------------------------------------------------------------------
    // YAML config (Phase-1): load once, validate lightly, and print summary
    // --------------------------------------------------------------------
    yamlcfg::Config cfg = yamlcfg::LoadConfig();
    if (const char* env = std::getenv("RJ_CLUSTER_UEPIPELINE"))
    {
        std::string s(env);
        if (s == "0" || s == "false" || s == "noSub") cfg.clusterUEpipeline = "noSub";
        else if (s == "1" || s == "true" || s == "variantA") cfg.clusterUEpipeline = "variantA";
        else if (s == "baseVariant") cfg.clusterUEpipeline = "baseVariant";
        else if (s == "variantB") cfg.clusterUEpipeline = "variantB";
        else cfg.clusterUEpipeline = s;  // pass through unknown for diagnostics
    }
    if (const char* env = std::getenv("RJ_PRESELECTION_VARIANT"))
    {
        std::string s = yamlcfg::NormalizePreselectionMode(std::string(env));
        if (!s.empty()) cfg.preselection = s;
    }
    if (const char* env = std::getenv("RJ_TIGHT_VARIANT"))
    {
        std::string s = yamlcfg::NormalizeTightMode(std::string(env));
        if (!s.empty()) cfg.tight = s;
    }
    if (const char* env = std::getenv("RJ_NONTIGHT_VARIANT"))
    {
        std::string s = yamlcfg::NormalizeNonTightMode(std::string(env));
        if (!s.empty()) cfg.nonTight = s;
    }
    if (!yamlcfg::IsPreselectionMode(cfg.preselection))
    {
        detail::bail("preselection must be 'reference', 'newPPG12', 'noPreCriteria', 'onlyNPB', 'refPlusNPB', or 'auauOnlyNPB'. Old variantA/B/C/D/E aliases are accepted in YAML/env parsing.");
    }
    if (!yamlcfg::IsTightMode(cfg.tight))
    {
        detail::bail("tight must be 'reference', 'newPPG12', 'auauEmbeddedBDT', 'centINDcontrol', 'centAsFeat', 'centDepBDTs', or one of the auau* validation BDT modes.");
    }
    if (!yamlcfg::IsNonTightMode(cfg.nonTight))
    {
        detail::bail("nonTight must be 'reference', 'newPPG12', 'auauBDTSideband', 'auauBDTComplement', 'centINDcontrol', 'centAsFeat', or 'centDepBDTs'.");
    }

    std::vector<idfanout::Entry> idFanoutEntries = idfanout::LoadFromEnv();
    if (idFanoutEntries.empty())
    {
        idfanout::Entry single;
        single.outRoot = outRoot;
        single.cfgTag = "";
        single.preselection = cfg.preselection;
        single.tight = cfg.tight;
        single.nonTight = cfg.nonTight;
        idFanoutEntries.push_back(single);
    }

    auto envFlag = [](const char* key) -> bool
    {
        const char* raw = std::getenv(key);
        if (!raw) return false;
        std::string v(raw);
        std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c){ return std::tolower(c); });
        return v == "1" || v == "true" || v == "yes" || v == "on";
    };
    const bool ppg12TableQAEnabled = envFlag("RJ_PPG12_TABLE_QA");
    const bool ppg12TableQAWantsPPScoreNodes =
        ppg12TableQAEnabled && !isAuAuRequested && !isSimEmbedded;
    const bool ppg12PhotonYieldPPSimWantsPPScoreNodes =
        envFlag("RJ_PPG12_PHOTON_YIELD") && isSim && !isAuAuRequested && !isSimEmbedded;

    bool fanoutUsesNPB = false;
    bool fanoutUsesAuAuNPB = false;
    bool fanoutUsesNewPPG12Tight = false;
    for (const auto& e : idFanoutEntries)
    {
        fanoutUsesNPB = fanoutUsesNPB || yamlcfg::PreselectionUsesNPB(e.preselection);
        fanoutUsesAuAuNPB = fanoutUsesAuAuNPB || yamlcfg::PreselectionUsesAuAuNPB(e.preselection);
        fanoutUsesNewPPG12Tight = fanoutUsesNewPPG12Tight || (e.tight == "newPPG12");
    }
    // PPG12 table-QA needs the analysis-row BDT score and the ordinary NPB
    // score. It should not force an AuAu-only NPB model unless the selected
    // photon-ID row itself uses preselection=auauOnlyNPB.
    fanoutUsesNewPPG12Tight =
        fanoutUsesNewPPG12Tight ||
        ppg12TableQAWantsPPScoreNodes ||
        ppg12PhotonYieldPPSimWantsPPScoreNodes;
    const bool ppPhotonIDTrainingWantsNPBAudit =
        isSim && !isSimEmbedded &&
        (cfg.pp_photonid_extract_only || cfg.pp_photonid_training_tree) &&
        !cfg.npb_model_file.empty() &&
        !cfg.npb_features.empty();
    const bool attachPPNPBScore =
        fanoutUsesNPB ||
        ppPhotonIDTrainingWantsNPBAudit ||
        ppg12TableQAWantsPPScoreNodes ||
        ppg12PhotonYieldPPSimWantsPPScoreNodes;

    cfg.preselection = idFanoutEntries.front().preselection;
    cfg.tight = idFanoutEntries.front().tight;
    cfg.nonTight = idFanoutEntries.front().nonTight;

    const std::vector<std::string> activeJetRKeys = yamlcfg::LoadJetRKeys(vlevel);
    
    std::vector<double> unfoldJetPtEdges;
    yamlcfg::ExpandUniformEdges(unfoldJetPtEdges,
                                cfg.unfold_jet_pt_start,
                                cfg.unfold_jet_pt_stop,
                                cfg.unfold_jet_pt_step);
    
    if (vlevel > 0)
    {
        std::cout << "\n[CFG] analysis_config.yaml\n"
        << "  path: " << cfg.yamlPath << "\n"
        << "  photon_eta_abs_max: " << cfg.photon_eta_abs_max << "\n"
        << "  jet_pt_min: " << cfg.jet_pt_min << "\n"
        << "  back_to_back_dphi_min_pi_fraction: " << cfg.back_to_back_dphi_min_pi_fraction
        << "  -> radians=" << (cfg.back_to_back_dphi_min_pi_fraction * M_PI) << "\n"
        << "  use_vz_cut: " << (cfg.use_vz_cut ? "true" : "false") << "\n"
        << "  vz_cut_cm: " << cfg.vz_cut_cm << "\n"
        << "  setMinBiasClassifer: " << (cfg.setMinBiasClassifer ? "true" : "false") << "\n"
        << "  coneR: " << cfg.isoConeR << "\n"
        << "  matching: {pho_dr_max=" << cfg.pho_dr_max << ", jet_dr_max=" << cfg.jet_dr_max << "}\n"
        << "  isolation_wp: {aGeV=" << cfg.isoA << ", bPerGeV=" << cfg.isoB
        << ", sideGapGeV=" << cfg.isoGap << ", fixedGeV=" << cfg.isoFixed
        << ", coneR=" << cfg.isoConeR
        << ", towerMin=" << cfg.isoTowMin
        << ", isSlidingIso=" << (cfg.isSlidingIso ? "true" : "false") << "}\n"
        << "  isSlidingAndFixed: " << (cfg.isSlidingAndFixed ? "true" : "false") << "\n"
        << "  fixedGeV: " << cfg.isoFixed << "\n"
        << "  auau_cent_iso_wp: " << cfg.auauCentIsoWP.size() << " entries\n";
        std::cout << "  pp_iso_wp_r30: " << (cfg.ppIsoWPR30.configured ? "configured" : "legacy fallback");
        if (cfg.ppIsoWPR30.configured)
        {
            std::cout << " {aGeV=" << cfg.ppIsoWPR30.aGeV
                      << ", bPerGeV=" << cfg.ppIsoWPR30.bPerGeV
                      << ", sideGapGeV=" << cfg.ppIsoWPR30.sideGapGeV << "}";
        }
        std::cout << "\n  pp_iso_wp_r40: " << (cfg.ppIsoWPR40.configured ? "configured" : "legacy fallback");
        if (cfg.ppIsoWPR40.configured)
        {
            std::cout << " {aGeV=" << cfg.ppIsoWPR40.aGeV
                      << ", bPerGeV=" << cfg.ppIsoWPR40.bPerGeV
                      << ", sideGapGeV=" << cfg.ppIsoWPR40.sideGapGeV << "}";
        }
        std::cout << "\n";
        std::cout << "  jes3_photon_pt_bins: [";
        for (std::size_t i = 0; i < cfg.jes3_photon_pt_bins.size(); ++i)
        {
            std::cout << cfg.jes3_photon_pt_bins[i] << (i + 1 < cfg.jes3_photon_pt_bins.size() ? ", " : "");
        }
        std::cout << "]\n  unfold_reco_photon_pt_bins: [";
        for (std::size_t i = 0; i < cfg.unfold_reco_photon_pt_bins.size(); ++i)
        {
            std::cout << cfg.unfold_reco_photon_pt_bins[i] << (i + 1 < cfg.unfold_reco_photon_pt_bins.size() ? ", " : "");
        }
        std::cout << "]\n  unfold_truth_photon_pt_bins: [";
        for (std::size_t i = 0; i < cfg.unfold_truth_photon_pt_bins.size(); ++i)
        {
            std::cout << cfg.unfold_truth_photon_pt_bins[i] << (i + 1 < cfg.unfold_truth_photon_pt_bins.size() ? ", " : "");
        }
        std::cout << "]\n  unfold_jet_pt_edges: start=" << cfg.unfold_jet_pt_start
        << " stop=" << cfg.unfold_jet_pt_stop
        << " step=" << cfg.unfold_jet_pt_step
        << " (nedges=" << unfoldJetPtEdges.size() << ")\n"
        << "  unfold_xj_bins: [";
        for (std::size_t i = 0; i < cfg.unfold_xj_bins.size(); ++i)
        {
            std::cout << cfg.unfold_xj_bins[i] << (i + 1 < cfg.unfold_xj_bins.size() ? ", " : "");
        }
        std::cout << "]\n"
        << "  leading_response_family: "
        << (cfg.leading_response_family.empty() ? "nominal" : cfg.leading_response_family)
        << "\n"
        << "  require_towerinfo_truth_matching: "
        << (cfg.require_towerinfo_truth_matching ? "true" : "false")
        << "\n"
        << "  clusterUEpipeline: " << cfg.clusterUEpipeline << "\n"
        << "  doPi0Analysis: " << (cfg.doPi0Analysis ? "true" : "false") << "\n"
        << "  event_display_tree: " << (cfg.event_display_tree ? "true" : "false") << "\n"
        << "  event_display_tree_max_per_bin: " << cfg.event_display_tree_max_per_bin << "\n\n";
        
        std::cout << "[CFG] caloInputMode: " << caloInputMode << "\n";
    }
    
    // --------------------------------------------------------------------
    // 2.  CDB + IO managers
    // --------------------------------------------------------------------
    recoConsts* rc = recoConsts::instance();
    const bool usePPG12PPSimRebuildCaloFromG4 =
        isSim && !isSimEmbedded &&
        env_truthy_local("RJ_PPG12_PHOTON_YIELD") &&
        env_truthy_local("RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4");
    const bool ppg12ClosureCanary =
        env_truthy_local("RJ_PPG12_CLOSURE_CANARY");
    const bool ppg12DINeutralityCanary =
        env_truthy_local("RJ_REPLAY_FOUNDATION_DI_NEUTRALITY_CANARY");
    const std::string ppg12ClosureCanaryId =
        std::getenv("RJ_PPG12_CLOSURE_CANARY_ID")
            ? detail::trim(std::string(std::getenv("RJ_PPG12_CLOSURE_CANARY_ID")))
            : std::string();
    const std::string ppg12DINeutralityCanaryId =
        std::getenv("RJ_REPLAY_FOUNDATION_DI_NEUTRALITY_CANARY_ID")
            ? detail::trim(std::string(std::getenv(
                  "RJ_REPLAY_FOUNDATION_DI_NEUTRALITY_CANARY_ID")))
            : std::string();
    const std::string ppg12HistoricalSeedSequence =
        "2991264730,4256268992,2394322166,874466025,2240380304";
    const int ppg12ExpectedPedestalSequence = 534;
    const std::string ppg12ReplaySeedSequence =
        std::getenv("RJ_PPG12_PPSIM_REPLAY_SEEDS")
            ? detail::trim(std::string(std::getenv("RJ_PPG12_PPSIM_REPLAY_SEEDS")))
            : std::string();
    const std::string ppg12ExpectedPedestalToken =
        std::getenv("RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE")
            ? detail::trim(std::string(
                  std::getenv("RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE")))
            : std::string();
    int ppg12NaturalPedestalSequence = -1;
    std::string ppg12PedestalFileForProvenance;

    if ((!ppg12ReplaySeedSequence.empty() ||
         !ppg12ExpectedPedestalToken.empty()) &&
        !ppg12ClosureCanary && !ppg12DINeutralityCanary)
    {
        detail::bail(
            "PPG12 historical RNG replay controls require the exact closure "
            "or replay-foundation DI-neutrality canary");
    }
    if (ppg12ClosureCanary && ppg12DINeutralityCanary)
        detail::bail("PPG12 closure and DI-neutrality canaries are mutually exclusive");
    if (ppg12DINeutralityCanary)
    {
        const bool safeCanaryId =
            !ppg12DINeutralityCanaryId.empty() &&
            ppg12DINeutralityCanaryId.size() <= 128 &&
            ppg12DINeutralityCanaryId.find_first_not_of(
                "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.:-") ==
                std::string::npos;
        if (!safeCanaryId)
            detail::bail("replay-foundation DI-neutrality canary ID is missing or unsafe");
        if (!env_truthy_local("RJ_REPLAY_FOUNDATION_CANARY") ||
            !isSim || isSimEmbedded || isAuAuRequested ||
            !usePPG12PPSimRebuildCaloFromG4 ||
            !env_truthy_local("RJ_PPG12_PPSIM_G4_ONLY") ||
            !env_truthy_local("RJ_PPG12_PHOTON_YIELD_DOUBLE") ||
            !env_truthy_local("RJ_PPG12_PERIOD_STRICT_DI"))
        {
            detail::bail(
                "replay-foundation DI-neutrality RNG controls require the "
                "pp-only archived double-interaction G4 rebuild canary path");
        }
        if (ppg12ReplaySeedSequence != ppg12HistoricalSeedSequence ||
            ppg12ExpectedPedestalToken !=
                std::to_string(ppg12ExpectedPedestalSequence))
        {
            detail::bail(
                "replay-foundation DI-neutrality canary requires the exact "
                "historical five-seed FIFO and pedestal sequence 534");
        }
        if (rc->FlagExist("RANDOMSEED"))
            detail::bail("DI-neutrality canary forbids recoConsts RANDOMSEED");
        PHRandomSeed::Verbosity(1);
        for (const unsigned int seed : {
                 2991264730U, 4256268992U, 2394322166U,
                 874466025U, 2240380304U})
            PHRandomSeed::LoadSeed(seed);
        std::cout << "[REPLAY_FOUNDATION_DI_NEUTRALITY_RNG] canary_id="
                  << ppg12DINeutralityCanaryId
                  << " mode=historical_fifo_replay_v2"
                  << " replay_sequence=" << ppg12ReplaySeedSequence
                  << " RANDOMSEED_absent=1" << std::endl;
    }
    if (ppg12ClosureCanary)
    {
        const bool safeCanaryId =
            !ppg12ClosureCanaryId.empty() &&
            ppg12ClosureCanaryId.size() <= 128 &&
            ppg12ClosureCanaryId.find_first_not_of(
                "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.:-") ==
                std::string::npos;
        if (!safeCanaryId)
            detail::bail("RJ_PPG12_CLOSURE_CANARY_ID is missing or unsafe");
        if (!isSim || isSimEmbedded || isAuAuRequested ||
            !usePPG12PPSimRebuildCaloFromG4 ||
            !env_truthy_local("RJ_PPG12_PPSIM_G4_ONLY"))
        {
            detail::bail(
                "PPG12 closure RNG controls require the pp-only PPG12 exact "
                "G4-only rebuild path");
        }
        if (std::getenv("RJ_PPG12_PPSIM_FIXED_RANDOMSEED") ||
            std::getenv("RJ_PPG12_PPSIM_FIXED_PEDESTAL_SEQUENCE"))
        {
            detail::bail(
                "synthetic fixed-seed controls are forbidden in the historical replay");
        }
        if (ppg12ReplaySeedSequence != ppg12HistoricalSeedSequence ||
            ppg12ExpectedPedestalToken !=
                std::to_string(ppg12ExpectedPedestalSequence))
        {
            detail::bail(
                "PPG12 closure canary requires the exact historical five-seed "
                "FIFO replay and pedestal sequence 534");
        }
        if (!env_truthy_local("RJ_DISABLE_JES_CDB_AUDIT"))
        {
            detail::bail(
                "PPG12 closure canary requires RJ_DISABLE_JES_CDB_AUDIT=1");
        }
        for (const char* forbidden : {
                 "RJ_CDB_GLOBALTAG",
                 "RJ_CDB_TIMESTAMP",
                 "RJ_TRUTH_JETS_MODE",
                 "RJ_DETAILED_CEMC_GEOM",
                 "RJ_CALO_INPUT_MODE",
                 "RJ_CALO_CLUSTER_NODE"})
        {
            if (std::getenv(forbidden))
            {
                detail::bail(
                    std::string("PPG12 closure canary forbids inherited override ") +
                    forbidden);
            }
        }
        if (std::getenv("RJ_PPG12_PEDESTAL_OVERRIDE") ||
            std::getenv("RJ_PPG12_DI_ARCHIVED_EXPECT_PEDESTAL"))
        {
            detail::bail(
                "generic or archived pedestal controls cannot be combined with the closure oracle");
        }
        if (rc->FlagExist("RANDOMSEED"))
        {
            detail::bail(
                "historical PPG12 seed replay requires recoConsts RANDOMSEED to be absent");
        }
        PHRandomSeed::Verbosity(1);
        std::cout << "[PPG12_CLOSURE_RNG] canary_id=" << ppg12ClosureCanaryId
                  << " mode=historical_fifo_replay_v2"
                  << " replay_sequence=" << ppg12ReplaySeedSequence
                  << " RANDOMSEED_absent=1" << std::endl;
    }

    // CDB_GLOBALTAG is REQUIRED for any CDBInterface::getUrl() call.
    // The PPG12 pp-SIM G4 rebuild macro uses MDC2; keep the normal analysis
    // default untouched outside that gated diagnostic/parity path.
    std::string gtag = usePPG12PPSimRebuildCaloFromG4 ? "MDC2" : "newcdbtag";
    if (!ppg12ClosureCanary)
    {
        if (const char* envGT = std::getenv("RJ_CDB_GLOBALTAG"))
        {
            std::string tmp = detail::trim(std::string(envGT));
            if (!tmp.empty()) gtag = tmp;
        }
    }
    rc->set_StringFlag("CDB_GLOBALTAG", gtag);
    
    if (vlevel > 0)
        std::cout << "[INFO] CDB_GLOBALTAG=" << gtag << "\n";

    // TIMESTAMP:
    //  - DATA: use run number (must be > 1000 so Calo_Calib treats it as DATA)
    //  - embedded SIM: use the embedded data run for CDB/centrality context.
    //    Blair-style DST_CALO embedded inputs still run Process_Calo_Calib() once
    //    to build TOWERINFO_CALIB_{CEMC,HCALIN,HCALOUT} and CLUSTERINFO_CEMC,
    //    while Calo_Calib skips data-only skimmer/ZDC/status pieces for embedded MC.
    //  - PPG12 pp-SIM G4 rebuild: use the SIM run number, matching PPG12's
    //    anatreemaker macro so Process_Calo_Calib() stays in its SIM branch.
    //  - plain SIM: use a known-good fixed timestamp (SIM does NOT run Process_Calo_Calib)
    unsigned long long cdbts = static_cast<unsigned long long>(run);
    
    if (!isSim || isSimEmbedded)
    {
        // DATA and embedded SIM carry real data-run context.
        if (cdbts <= 1000ULL)
        {
            std::cerr << "[FATAL] DATA/embedded run number " << cdbts
            << " is invalid for data-run CDB/centrality context.\n";
            throw std::runtime_error("Invalid DATA/embedded TIMESTAMP (must be > 1000).");
        }
    }
    else
    {
        cdbts = usePPG12PPSimRebuildCaloFromG4
            ? (run > 0 ? static_cast<unsigned long long>(run) : 28ULL)
            : 47289ULL;  // keep your old working SIM timestamp
    }
    
    if (!ppg12ClosureCanary)
    {
        if (const char* ts = std::getenv("RJ_CDB_TIMESTAMP"))
        {
            char* end = nullptr;
            unsigned long long tmp = std::strtoull(ts, &end, 10);
            if (end != ts && tmp > 0ULL) cdbts = tmp;
        }
    }
    else if (run != 28 || cdbts != 28ULL)
    {
        detail::bail(
            "PPG12 closure canary requires Run-28 input identity and TIMESTAMP=28");
    }
    
    rc->set_uint64Flag("TIMESTAMP", cdbts);
    
    if (vlevel > 0)
        std::cout << "[INFO] CDB TIMESTAMP=" << rc->get_uint64Flag("TIMESTAMP")
        << " (isSim=" << (isSim ? "true" : "false") << ")\n";
    
    CDBInterface::instance()->Verbosity(0);

    const bool auditJESCdb = !env_truthy_local("RJ_DISABLE_JES_CDB_AUDIT");
    if (auditJESCdb)
    {
        try
        {
            const std::string jesUrl = CDBInterface::instance()->getUrl("JES_Calib_Default");
            std::cout << "[INFO] JES_CDB_AUDIT CDB_GLOBALTAG=" << gtag
                      << " TIMESTAMP=" << rc->get_uint64Flag("TIMESTAMP")
                      << " JES_Calib_Default="
                      << (jesUrl.empty() ? std::string("<empty>") : jesUrl)
                      << "\n";
        }
        catch (const std::exception& e)
        {
            std::cout << "[WARN] JES_CDB_AUDIT failed for CDB_GLOBALTAG=" << gtag
                      << " TIMESTAMP=" << rc->get_uint64Flag("TIMESTAMP")
                      << ": " << e.what() << "\n";
        }
        catch (...)
        {
            std::cout << "[WARN] JES_CDB_AUDIT failed for CDB_GLOBALTAG=" << gtag
                      << " TIMESTAMP=" << rc->get_uint64Flag("TIMESTAMP")
                      << " with unknown exception\n";
        }
    }
    
    
    // The preserved PPG12 executable registers FlagHandler immediately after
    // InputInit/InputRegister.  Defer it only for the exact closure canary;
    // every ordinary pp/AuAu path keeps the established registration point.
    if (!ppg12ClosureCanary)
    {
        auto* flag = new FlagHandler();
        se->registerSubsystem(flag);
    }
    
    // ------------------------------------------------------------------
    // Decide how to source truth jets in SIM:
    //
    //   RJ_TRUTH_JETS_MODE=AUTO  (default)
    //       - if list has 3rd column (DST_JETS): read jets directly from DST
    //       - otherwise: build truth jets from TRUTH particles (TruthJetInput)
    //
    //   RJ_TRUTH_JETS_MODE=DST
    //       - require 3rd column (DST_JETS) and read truth jets from DST
    //
    //   RJ_TRUTH_JETS_MODE=BUILD
    //       - ignore 3rd column and build truth jets from TRUTH particles
    //
    //   RJ_TRUTH_JETS_MODE=BOTH
    //       - read DST jets AND also build a second "from particles" truth-jet
    //         collection under a different node name for QA comparisons
    // ------------------------------------------------------------------
    bool useDSTTruthJets = false;
    bool buildTruthJetsFromParticles = false;
    bool buildTruthJetsAsAltNode = false;   // only true in BOTH mode
    
    auto all_nonempty = [](const std::vector<std::string>& v) -> bool
    {
        if (v.empty()) return false;
        for (const auto& s : v) if (s.empty()) return false;
        return true;
    };
    
    std::map<std::string, std::string> frogResolvedCache;
    std::size_t nFrogResolved = 0;
    auto resolveInputPath = [&](const std::string& path, const char* streamLabel) -> std::string
    {
        if (path.empty() || path == "NONE") return path;
        if (path.find("://") != std::string::npos) return path;
        if (!path.empty() && path[0] == '/') return path;
        if (gSystem && !gSystem->AccessPathName(path.c_str())) return path;

        auto cached = frogResolvedCache.find(path);
        if (cached != frogResolvedCache.end()) return cached->second;

        FROG frog;
        const char* resolvedRaw = frog.location(path);
        if (resolvedRaw && std::string(resolvedRaw).size())
        {
            std::string resolved(resolvedRaw);
            if (resolved.find("://") == std::string::npos &&
                (resolved.empty() || resolved[0] != '/') &&
                (!gSystem || gSystem->AccessPathName(resolved.c_str())))
            {
                // Fun4AllDstInputManager has its own logical-file handling; do
                // not reject unresolved logical names here solely because this
                // direct FROG probe returned the original basename.
                frogResolvedCache[path] = path;
                if (vlevel > 1)
                {
                    std::cout << "[FROG] " << streamLabel
                              << ": leaving logical input for Fun4All: "
                              << path << std::endl;
                }
                return path;
            }
            frogResolvedCache[path] = resolved;
            if (resolved != path) ++nFrogResolved;
            if (vlevel > 1)
            {
                std::cout << "[FROG] " << streamLabel << ": " << path
                          << " -> " << resolved << std::endl;
            }
            return resolved;
        }

        frogResolvedCache[path] = path;
        if (vlevel > 1)
        {
            std::cout << "[FROG] " << streamLabel
                      << ": no direct resolution; leaving logical input for Fun4All: "
                      << path << std::endl;
        }
        return path;
    };

    auto resolveVectorPaths = [&](std::vector<std::string>& paths, const char* streamLabel)
    {
        for (auto& path : paths)
        {
            path = resolveInputPath(path, streamLabel);
        }
    };

    resolveVectorPaths(filesCalo, "CALO");
    resolveVectorPaths(filesZdc, "ZDC_RAW");
    resolveVectorPaths(filesG4, "G4");
    resolveVectorPaths(filesJets, "JETS");
    resolveVectorPaths(filesGlobal, "GLOBAL");
    resolveVectorPaths(filesMbd, "MBD_EPD");

    if ((verbose || vlevel > 0) && nFrogResolved > 0)
    {
        std::cout << "[FROG] resolved " << nFrogResolved
                  << " logical input file names before Fun4All AddFile" << std::endl;
    }

    const bool listHasZdc    = all_nonempty(filesZdc);
    const bool listHasCalo   = all_nonempty(filesCalo);
    const bool listHasG4     = all_nonempty(filesG4);
    const bool listHasJets   = all_nonempty(filesJets);
    const bool listHasGlobal = all_nonempty(filesGlobal);
    const bool listHasMbd    = all_nonempty(filesMbd);
    const std::string simSampleLower = env_lower("RJ_SIM_SAMPLE");
    const bool candidateSkimNeedsPhotonStitchTruth =
        auauCandidateSkimOnly && isSimEmbedded &&
        simSampleLower.find("embeddedphoton") != std::string::npos;
    const bool candidateSkimNeedsInclusiveStitchTruth =
        auauCandidateSkimOnly && isSimEmbedded &&
        !candidateSkimNeedsPhotonStitchTruth &&
        simSampleLower.find("embeddedjet") != std::string::npos;
    if (auauCandidateSkimOnly && isSimEmbedded &&
        !candidateSkimNeedsPhotonStitchTruth &&
        !candidateSkimNeedsInclusiveStitchTruth)
    {
        detail::bail(
            "RJ_AUAU_CANDIDATE_SKIM_ONLY=1 on embedded simulation requires "
            "RJ_SIM_SAMPLE to identify an embeddedPhoton or embeddedJet sample; "
            "the generator-slice ownership gate must not be silently bypassed.");
    }
    const bool isRun24PPData = !isSim && !isPPrun25 && !isAuAuRequested;
    const bool usePPG12PPDataPair =
        isRun24PPData &&
        listHasCalo &&
        listHasZdc &&
        env_bool_local("RJ_PPG12_PP_DATA_PAIRED", true);
    auto first_nonempty_lower = [](const std::vector<std::string>& paths) -> std::string
    {
        for (const auto& path : paths)
        {
            if (path.empty() || path == "NONE") continue;
            std::string out = path;
            std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c){ return std::tolower(c); });
            return out;
        }
        return {};
    };
    const std::string firstAuxLower = first_nonempty_lower(filesZdc);
    const bool auxLooksZdcRaw =
        firstAuxLower.find("zdc") != std::string::npos;
    const bool auxLooksJetCalo =
        firstAuxLower.find("jetcalo") != std::string::npos ||
        firstAuxLower.find("dst_jetcalo") != std::string::npos;
    const bool useAuAuJetCaloDataPair =
        isAuAuRequested &&
        !isSim &&
        listHasCalo &&
        listHasZdc &&
        env_bool_local("RJ_AUAU_DATA_PAIRED", true) &&
        (caloInputMode == "jetcalo" || auxLooksJetCalo) &&
        !auxLooksZdcRaw;
    const bool usePPG12PPSimG4OnlyInput =
        usePPG12PPSimRebuildCaloFromG4 &&
        env_truthy_local("RJ_PPG12_PPSIM_G4_ONLY");
    // Production-gated reconstruction arm for the deployed PPG12
    // double-interaction contract.  Keep the current RecoilJets binary and
    // release ABI, but reproduce the proven archived subsystem order:
    // G4Hits + DST truth jets, MBD/vertex reconstruction, standard tower
    // helpers, and Process_Calo_Calib with CaloTowerStatus enabled.  The
    // current run-28 helper publishes the raw simulated calorimeter containers
    // as TOWERS_*; Process_Calo_Calib then copies their calibrated/status state
    // into TOWERINFO_CALIB_*, matching the proven current-ABI canary.
    const bool requestPPG12ArchivedDIG4OnlyReco =
        env_truthy_local("RJ_PPG12_DI_ARCHIVED_RECO_CHAIN");
    if (requestPPG12ArchivedDIG4OnlyReco && !usePPG12PPSimG4OnlyInput)
    {
        detail::bail(
            "RJ_PPG12_DI_ARCHIVED_RECO_CHAIN=1 requires both "
            "RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1 and "
            "RJ_PPG12_PPSIM_G4_ONLY=1");
    }
    const bool usePPG12ArchivedDIG4OnlyReco =
        requestPPG12ArchivedDIG4OnlyReco && usePPG12PPSimG4OnlyInput;

    if (usePPG12ArchivedDIG4OnlyReco)
    {
        rc->set_IntFlag("RUNNUMBER", run);
        const std::string simSample = env_lower("RJ_SIM_SAMPLE");
        const std::string embeddedSample =
            env_lower("RJ_EMBEDDED_INCLUSIVE_JET_SAMPLE");
        const std::string archivedTruthMode =
            env_lower("RJ_TRUTH_JETS_MODE", "auto");
        const std::vector<std::string> allowedPhotonSamples = {
            "run28_photonjet5_double", "run28_photonjet10_double",
            "run28_photonjet20_double"};
        const std::vector<std::string> allowedInclusiveSamples = {
            "run28_jet8_double", "run28_jet12_double",
            "run28_jet20_double", "run28_jet30_double",
            "run28_jet40_double"};
        const bool run28PhotonDoubleSample =
            std::find(allowedPhotonSamples.begin(),
                      allowedPhotonSamples.end(), simSample) !=
            allowedPhotonSamples.end();
        const bool run28InclusiveDoubleSample =
            std::find(allowedInclusiveSamples.begin(),
                      allowedInclusiveSamples.end(), simSample) !=
            allowedInclusiveSamples.end();
        const bool run28DoubleSample =
            run28PhotonDoubleSample || run28InclusiveDoubleSample;
        const std::string sampleSlice = run28DoubleSample
            ? simSample.substr(std::string("run28_").size(),
                               simSample.size() - std::string("run28_").size() -
                                   std::string("_double").size())
            : std::string();
        auto all_paths_contain = [](const std::vector<std::string>& paths,
                                    const std::string& token) -> bool
        {
            if (paths.empty()) return false;
            for (const auto& path : paths)
            {
                if (path.empty() || path.find(token) == std::string::npos)
                    return false;
            }
            return true;
        };

        if ((run28PhotonDoubleSample && datasetToken != "issim") ||
            (run28InclusiveDoubleSample && datasetToken != "issiminclusive"))
        {
            detail::bail(
                "RJ_PPG12_DI_ARCHIVED_RECO_CHAIN=1 requires photonjet*_double "
                "with RJ_DATASET=isSim or jet*_double with "
                "RJ_DATASET=isSimInclusive");
        }
        if (!run28DoubleSample || embeddedSample != simSample)
        {
            detail::bail(
                "RJ_PPG12_DI_ARCHIVED_RECO_CHAIN=1 requires matching "
                "run28_photonjet{5,10,20}_double or "
                "run28_jet{8,12,20,30,40}_double sample identities");
        }
        if (archivedTruthMode != "dst")
        {
            detail::bail(
                "RJ_PPG12_DI_ARCHIVED_RECO_CHAIN=1 requires "
                "RJ_TRUTH_JETS_MODE=DST");
        }
        if (!env_truthy_local("RJ_PPG12_PHOTON_YIELD_DOUBLE") ||
            !env_truthy_local("RJ_PPG12_PERIOD_STRICT_DI") ||
            env_truthy_local("RJ_PPG12_PERIOD_ALLOW_ALL_SIM"))
        {
            detail::bail(
                "RJ_PPG12_DI_ARCHIVED_RECO_CHAIN=1 requires the strict "
                "double-interaction period contract");
        }
        if (run != 28 || rc->get_IntFlag("RUNNUMBER") != 28 ||
            rc->get_uint64Flag("TIMESTAMP") != 28ULL || gtag != "MDC2")
        {
            detail::bail(
                "RJ_PPG12_DI_ARCHIVED_RECO_CHAIN=1 requires the run-28 MDC2 "
                "CDB contract");
        }
        if (listHasCalo || listHasGlobal || listHasMbd ||
            !listHasG4 || !listHasJets)
        {
            detail::bail(
                "RJ_PPG12_DI_ARCHIVED_RECO_CHAIN=1 requires exactly "
                "NONE,G4Hits,DST_JETS,NONE,NONE input rows");
        }
        if (!all_paths_contain(
                filesG4, "/js_pp200_signal_dual/g4hits/run0028/" +
                             sampleSlice + "/") ||
            !all_paths_contain(
                filesJets, "/js_pp200_signal_dual/nopileup/jets/run0028/" +
                               sampleSlice + "/"))
        {
            detail::bail(
                "RJ_PPG12_DI_ARCHIVED_RECO_CHAIN=1 requires matching "
                "run-28 dual-interaction G4Hits and truth-jet sources");
        }
    }
    const bool usePPG12Fig11G4OnlyRebuild =
        usePPG12PPSimG4OnlyInput &&
        env_truthy_local("RJ_PPG12_FIG11_SB_DIAGNOSTIC") &&
        !usePPG12ArchivedDIG4OnlyReco;

    if (usePPG12PPSimRebuildCaloFromG4)
    {
        if (!listHasG4 || !listHasJets)
        {
            detail::bail(
                "RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1 requires SIM list columns "
                "<DST_CALO_CLUSTER_or_NONE> <G4Hits> <DST_JETS> [<DST_GLOBAL_or_NONE> <DST_MBD_EPD_or_NONE>], "
                "matching the PPG12 anatreemaker input contract.");
        }
        if (!usePPG12PPSimG4OnlyInput && (!listHasCalo || !listHasMbd))
        {
            detail::bail(
                "RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1 requires the PPG12 SIM input stack "
                "lanes 1=DST_CALO_CLUSTER and 3=DST_MBD_EPD in addition to lanes 0=G4Hits and 4=DST_JETS. "
                "Set RJ_PPG12_PPSIM_G4_ONLY=1 for the live executable PPG12 "
                "SI/DI contract that intentionally provides only G4Hits + DST_JETS.");
        }
        if (ppg12ClosureCanary &&
            (!usePPG12PPSimG4OnlyInput || listHasCalo || listHasGlobal || listHasMbd))
        {
            detail::bail(
                "The executable PPG12 closure oracle requires exactly G4Hits + "
                "DST_TRUTH_JET inputs; CALO, GLOBAL, and MBD lanes must be NONE");
        }
    }

    const bool needZdcRawForMinBias =
        cfg.setMinBiasClassifer &&
        !useAuAuJetCaloDataPair &&
        !isSim &&
        (isAuAuRequested || (!isPPrun25 && run > 53864));

    if (needZdcRawForMinBias && !listHasZdc)
    {
        detail::bail("setMinBiasClassifer=true for AuAu data requires a paired DST_ZDC_RAW file in column 2 of each input-list row.");
    }
    
    if (verbose || vlevel > 0)
    {
        std::cout << "[FLOW] input contract:"
        << " caloFiles=" << filesCalo.size()
        << " | ZDC_RAW=" << (listHasZdc ? "present" : "missing")
        << " | G4=" << (listHasG4 ? "present" : "missing")
        << " | DST_JETS=" << (listHasJets ? "present" : "missing")
        << " | DST_GLOBAL=" << (listHasGlobal ? "present" : "missing")
        << " | DST_MBD_EPD=" << (listHasMbd ? "present" : "missing")
        << " | caloInputMode=" << caloInputMode
        << " | ppPPG12DataPair=" << (usePPG12PPDataPair ? "true" : "false")
        << " | auauJetCaloDataPair=" << (useAuAuJetCaloDataPair ? "true" : "false")
        << " | minBiasClassifierGate=" << (cfg.setMinBiasClassifer ? "true" : "false")
        << std::endl;
    }
    
    const std::string truthMode = env_lower("RJ_TRUTH_JETS_MODE", "auto");
    if (truthMode == "dst")
    {
        useDSTTruthJets = true;
        buildTruthJetsFromParticles = false;
        buildTruthJetsAsAltNode = false;
    }
    else if (truthMode == "build")
    {
        useDSTTruthJets = false;
        buildTruthJetsFromParticles = true;
        buildTruthJetsAsAltNode = false;
    }
    else if (truthMode == "both")
    {
        useDSTTruthJets = true;
        buildTruthJetsFromParticles = true;
        buildTruthJetsAsAltNode = true;
    }
    else
    {
        // AUTO (or any unrecognized token): prefer DST truth jets if provided
        useDSTTruthJets = listHasJets;
        buildTruthJetsFromParticles = !listHasJets;
        buildTruthJetsAsAltNode = false;
    }

    if (auauCandidateSkimOnly)
    {
        // The bounded shower-contract diagnostic consumes reconstructed
        // calorimeter candidates, but embedded sample ownership is still a
        // physics gate inside RecoilJets.  Retain only the truth lane needed
        // by that gate: G4/HepMC for Photon12/20, or DST truth jets for
        // Jet12/20/30/40.  Do not rebuild any alternate truth-jet collection.
        // This remains strictly environment-gated and leaves production
        // steering unchanged.
        useDSTTruthJets = candidateSkimNeedsInclusiveStitchTruth;
        buildTruthJetsFromParticles = false;
        buildTruthJetsAsAltNode = false;

        if (candidateSkimNeedsPhotonStitchTruth && !listHasG4)
        {
            detail::bail(
                "Embedded-photon candidate skim requires the G4/HepMC input "
                "lane for the Photon12/20 generator-slice ownership gate.");
        }
        if (candidateSkimNeedsInclusiveStitchTruth && !listHasJets)
        {
            detail::bail(
                "Embedded-inclusive candidate skim requires the DST_JETS input "
                "lane for the Jet12/20/30/40 truth-jet ownership gate.");
        }
    }
    
    if (usePPG12PPSimRebuildCaloFromG4)
    {
#if RJ_HAS_SPHENIX_G4_INPUT_MACROS
        auto write_ppg12_stack_list =
            [&](const std::vector<std::string>& paths,
                const std::string& label) -> std::string
        {
            const char* dirEnv = std::getenv("RJ_PPG12_PPSIM_INPUT_STACK_DIR");
            std::string dir = (dirEnv && std::string(dirEnv).size())
                ? detail::trim(std::string(dirEnv))
                : std::string("/tmp");
            if (dir.empty()) dir = "/tmp";
            if (!dir.empty() && dir.back() == '/') dir.pop_back();

            std::ostringstream path;
            path << dir << "/rj_ppg12_" << label << "_"
                 << static_cast<long long>(::getpid()) << ".list";
            std::ofstream out(path.str());
            if (!out.is_open())
            {
                detail::bail("failed to create PPG12 SIM input-stack list: " + path.str());
            }
            for (const auto& p : paths)
            {
                if (!p.empty()) out << p << '\n';
            }
            out.close();
            return path.str();
        };

        const bool usePPG12PPSimAuxInputs = !usePPG12PPSimG4OnlyInput;
        const std::string g4List   = write_ppg12_stack_list(filesG4, "g4hits");
        const std::string caloList = (usePPG12PPSimAuxInputs && listHasCalo)
            ? write_ppg12_stack_list(filesCalo, "dst_calo_cluster")
            : std::string();
        const std::string mbdList = (usePPG12PPSimAuxInputs && listHasMbd)
            ? write_ppg12_stack_list(filesMbd, "dst_mbd_epd")
            : std::string();
        const std::string jetsList = write_ppg12_stack_list(filesJets, "dst_truth_jet");

        Input::VERBOSITY = (vlevel > 0) ? 1 : 0;
        Input::READHITS = true;
        INPUTREADHITS::listfile[0] = g4List;
        if (usePPG12PPSimAuxInputs && listHasCalo) INPUTREADHITS::listfile[1] = caloList;
        if (usePPG12PPSimAuxInputs && listHasMbd) INPUTREADHITS::listfile[3] = mbdList;
        INPUTREADHITS::listfile[4] = jetsList;
        if (ppg12ClosureCanary)
        {
            // The executable oracle includes the Calo_Calib status chain.
            unsetenv("RJ_SKIP_CALO_TOWER_STATUS");
            unsetenv("RJ_CALO_TOWER_STATUS_INPUT_PREFIX");
        }
        else if (usePPG12PPSimG4OnlyInput &&
                 !usePPG12ArchivedDIG4OnlyReco)
        {
            // Preserve the pre-existing non-canary G4-rebuild behavior.
            setenv("RJ_SKIP_CALO_TOWER_STATUS", "1", 1);
        }
        else if (usePPG12ArchivedDIG4OnlyReco)
        {
            unsetenv("RJ_SKIP_CALO_TOWER_STATUS");
            // CEMC_Towers/HCAL*_Towers publish TOWERS_* in this G4 waveform
            // graph.  Fail closed against a stale external prefix override;
            // CaloTowerStatus and CaloTowerCalib must consume those nodes.
            unsetenv("RJ_CALO_TOWER_STATUS_INPUT_PREFIX");
        }
        else
        {
            unsetenv("RJ_SKIP_CALO_TOWER_STATUS");
        }
        InputInit();
        InputRegister();

        if (ppg12ClosureCanary)
        {
            // Exact live PPG12 graph for both SI and DI:
            // InputRegister -> FlagHandler -> Mbd_Reco helper ->
            // GlobalVertexReco -> RunSettings/towers -> InputManagers.
            auto* ppg12Flag = new FlagHandler();
            se->registerSubsystem(ppg12Flag);

            Enable::MBDRECO = true;
            Mbd_Reco();

            if (vlevel > 0)
            {
                std::cout << "[PPG12_PPSIM_REBUILD_CALO_FROM_G4] registering "
                          << "Mbd_Reco() + GlobalVertexReco after InputRegister() "
                          << "to match the live PPG12 SI/DI executable graph"
                          << std::endl;
            }
            std::unique_ptr<GlobalVertexReco> gvertex = std::make_unique<GlobalVertexReco>();
            gvertex->Verbosity(0);
            se->registerSubsystem(gvertex.release());

            RunSettings(28);
            Enable::CEMC_TOWERINFO = true;
            Enable::HCALIN_TOWERINFO = true;
            Enable::HCALOUT_TOWERINFO = true;
            CEMC_Towers();
            HCALInner_Towers();
            HCALOuter_Towers();

            auto* timerStats = new TimerStats();
            timerStats->OutFileName("jobtime.root");
            se->registerSubsystem(timerStats);
        }
        else
        {
            Enable::MBDRECO = false;
            if (usePPG12PPSimG4OnlyInput)
            {
                if (vlevel > 0)
                {
                    std::cout << "[PPG12_PPSIM_REBUILD_CALO_FROM_G4] registering early "
                              << "MbdDigitization + MbdReco + GlobalVertexReco "
                              << "for the non-canary G4-only rebuild" << std::endl;
                }
                std::unique_ptr<MbdDigitization> mbddigi =
                    std::make_unique<MbdDigitization>();
                se->registerSubsystem(mbddigi.release());
            }

            if (!usePPG12PPSimG4OnlyInput && vlevel > 0)
            {
                std::cout << "[PPG12_PPSIM_REBUILD_CALO_FROM_G4] registering "
                          << "MbdReco + GlobalVertexReco before Process_Calo_Calib "
                          << "for four-lane Fig.11-compatible input" << std::endl;
            }

            std::unique_ptr<MbdReco> mbdreco = std::make_unique<MbdReco>();
            se->registerSubsystem(mbdreco.release());

            std::unique_ptr<GlobalVertexReco> gvertex = std::make_unique<GlobalVertexReco>();
            se->registerSubsystem(gvertex.release());

            if (usePPG12PPSimG4OnlyInput)
            {
                RunSettings(28);
                Enable::CEMC_TOWERINFO = true;
                Enable::HCALIN_TOWERINFO = true;
                Enable::HCALOUT_TOWERINFO = true;
                if (usePPG12Fig11G4OnlyRebuild)
                {
                    CEMC_Cells();
                    HCALInner_Cells();
                    HCALOuter_Cells();
                }
                CEMC_Towers();
                HCALInner_Towers();
                HCALOuter_Towers();

                if (usePPG12Fig11G4OnlyRebuild)
                {
                    auto* clusterBuilder =
                        new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
                    clusterBuilder->Detector("CEMC");
                    clusterBuilder->set_threshold_energy(0.070);
                    const char* calibrationRoot = std::getenv("CALIBRATIONROOT");
                    if (!calibrationRoot || !std::string(calibrationRoot).size())
                        detail::bail("CALIBRATIONROOT is required for the PPG12 cluster builder");
                    clusterBuilder->LoadProfile(
                        std::string(calibrationRoot) +
                        "/EmcProfile/CEMCprof_Thresh30MeV.root");
                    clusterBuilder->set_UseTowerInfo(1);
                    se->registerSubsystem(clusterBuilder);
                }
            }
        }
        if (!usePPG12PPSimG4OnlyInput && (verbose || vlevel > 0))
        {
            std::cout << "[PPG12_FIG11_SB][FOUR_LANE] using DST_CALO_CLUSTER + DST_MBD_EPD inputs; "
                      << "not registering the G4 waveform/tower helper stack. Process_Calo_Calib() will match the visible PPG12 macro."
                      << std::endl;
            }
        InputManagers();

        TRandom3 randGen;
        randGen.SetSeed(PHRandomSeed());
        // Preserve the executable PPG12 random-call graph.  In the exact
        // closure oracle an external wrapper preloads the five seeds captured
        // in historical OutDir0, so the natural result must itself be 00534;
        // no pedestal identity is substituted.
        const int naturalSequence = randGen.Integer(3260);
        const int sequence = naturalSequence;
        ppg12NaturalPedestalSequence = naturalSequence;
        if ((ppg12ClosureCanary || ppg12DINeutralityCanary) &&
            naturalSequence != ppg12ExpectedPedestalSequence)
        {
            detail::bail(
                "historical PPG12 seed replay did not reproduce pedestal sequence 534");
        }
        std::ostringstream pedName;
        pedName << "pedestal-54256-0" << std::setw(4) << std::setfill('0') << sequence << ".root";
        ppg12PedestalFileForProvenance = pedName.str();
        auto* pedIn = new Fun4AllNoSyncDstInputManager("DST2");
        pedIn->AddFile(pedName.str());
        pedIn->Repeat();
        se->registerInputManager(pedIn);
        permittedRepeatingPedestalInputManager = pedIn;

        if (ppg12ClosureCanary)
        {
            std::cout << "[PPG12_CLOSURE_PEDESTAL] canary_id="
                      << ppg12ClosureCanaryId
                      << " natural_sequence=" << naturalSequence
                      << " chosen_sequence=" << sequence
                      << " logical_file=" << pedName.str()
                      << " phrandomseed_call_consumed=1" << std::endl;
        }

        if (verbose || vlevel > 0)
        {
            std::cout << "[PPG12_PPSIM_REBUILD_CALO_FROM_G4] using PPG12-style InputRegister/InputManagers"
                      << " g4List=" << g4List
                      << " caloList=" << ((usePPG12PPSimAuxInputs && listHasCalo) ? caloList : "DISABLED")
                      << " mbdList=" << ((usePPG12PPSimAuxInputs && listHasMbd) ? mbdList : "DISABLED")
                      << " jetsList=" << jetsList
                      << " listfile_indices=0"
                      << ((usePPG12PPSimAuxInputs && listHasCalo) ? ",1" : "")
                      << ((usePPG12PPSimAuxInputs && listHasMbd) ? ",3" : "")
                      << ",4"
                      << " g4_only=" << (usePPG12PPSimG4OnlyInput ? "true" : "false")
                      << " pedestal=" << pedName.str()
                      << std::endl;
        }
#else
        detail::bail(
            "RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1 requested, but GlobalVariables.C/G4_Input.C "
            "were not available to this ROOT macro environment.");
#endif
    }
    else
    {
        const bool usePPG12PPSimPrebuiltG4OnlyInput =
            isSim && !isSimEmbedded &&
            env_truthy_local("RJ_PPG12_FIG11_SB_DIAGNOSTIC") &&
            env_truthy_local("RJ_PPG12_PPSIM_G4_ONLY") &&
            env_truthy_local("RJ_ALLOW_LEGACY_PPG12_PPSIM_PREBUILT_G4_ONLY");

        if (isSim && !isSimEmbedded &&
            env_truthy_local("RJ_PPG12_FIG11_SB_DIAGNOSTIC") &&
            env_truthy_local("RJ_PPG12_PPSIM_G4_ONLY") &&
            !usePPG12PPSimPrebuiltG4OnlyInput)
        {
            detail::bail(
                "RJ_PPG12_FIG11_SB_DIAGNOSTIC=1 with RJ_PPG12_PPSIM_G4_ONLY=1 "
                "must use the PPG12-style G4 rebuild stack. Set "
                "RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1 so the macro runs "
                "Input::READHITS/InputRegister, Mbd_Reco, GlobalVertexReco, "
                "and the calo reconstruction path instead of the legacy "
                "prebuilt-G4 shortcut. To debug the old shortcut explicitly, "
                "set RJ_ALLOW_LEGACY_PPG12_PPSIM_PREBUILT_G4_ONLY=1.");
        }

        if (usePPG12PPSimPrebuiltG4OnlyInput)
        {
            // Legacy diagnostic-only path.  PPG12 production macros rebuild
            // from G4Hits with Input::READHITS/InputRegister, Mbd_Reco,
            // GlobalVertexReco, and the calo reconstruction chain; keep this
            // shortcut opt-in so it cannot become the canonical parity path.
            if (!listHasG4 || !listHasJets)
            {
                detail::bail(
                    "RJ_PPG12_FIG11_SB_DIAGNOSTIC=1 with RJ_PPG12_PPSIM_G4_ONLY=1 "
                    "requires G4Hits and DST_JETS columns for the prebuilt-cluster input mode.");
            }

            auto* inG4 = new Fun4AllNoSyncDstInputManager("DST_G4HITS_IN");
            for (const auto& f : filesG4) inG4->AddFile(f);
            se->registerInputManager(inG4);

            if (useDSTTruthJets)
            {
                auto* inJets = new Fun4AllNoSyncDstInputManager("DST_JETS_IN");
                for (const auto& f : filesJets) inJets->AddFile(f);
                se->registerInputManager(inJets);
            }

            if (vlevel > 0)
            {
                std::cout << "[PPG12_FIG11_SB][G4_ONLY_PREBUILT] registered G4Hits"
                          << (useDSTTruthJets ? " + DST_JETS" : "")
                          << " no-sync input streams without DST_GLOBAL and without calo rebuild"
                          << " (nFiles=" << filesG4.size() << ")" << std::endl;
            }
        }
        else
        {
        if (usePPG12PPDataPair || useAuAuJetCaloDataPair)
        {
            auto* inDataJet = new Fun4AllDstInputManager("DST_JET_IN");
            // CentralityReco must create its own v2 node when recalculating
            // AuAu DATA. Importing the old DST's v1 node shadows that node and
            // silently loses the new integer-bin witness (v1 has no bin API).
            // Only the recomputed object is suppressed; MBD inputs are kept.
            const bool replaceInputCentrality = useAuAuJetCaloDataPair &&
                !rj_centrality_binding::env("RJ_AUAU_CENTRALITY_SOURCE").empty();
            if (replaceInputCentrality)
                inDataJet->BranchSelect("*CentralityInfo", 0);
            for (const auto& f : filesCalo) inDataJet->AddFile(f);
            se->registerInputManager(inDataJet);

            auto* inDataJetCalo = new Fun4AllDstInputManager("DST_JETCALO_IN");
            if (replaceInputCentrality)
                inDataJetCalo->BranchSelect("*CentralityInfo", 0);
            for (const auto& f : filesZdc) inDataJetCalo->AddFile(f);
            se->registerInputManager(inDataJetCalo);
            if (replaceInputCentrality)
                std::cout << "CENTRALITY_INPUT action=RECOMPUTE_V2"
                          << " suppressed=CentralityInfo inputs=DST_JET,DST_JETCALO"
                          << " mbd_inputs=PRESERVED" << std::endl;

            if (vlevel > 0)
            {
                std::cout << (useAuAuJetCaloDataPair
                              ? "[INFO] AuAu data paired input enabled: registered DST_JET + DST_JETCALO streams"
                              : "[INFO] PPG12 pp-data paired input enabled: registered DST_Jet + DST_JETCALO streams")
                          << " (nFiles=" << filesCalo.size() << ")" << std::endl;
            }
        }
        else
        {
            // ------------------ Calo cluster DST (single-stream data or SIM calo lane) -------------------
            auto* inCalo = new Fun4AllDstInputManager("DSTcalofitting");
            for (const auto& f : filesCalo) inCalo->AddFile(f);
            se->registerInputManager(inCalo);

            // ------------------ ZDC RAW DST (AuAu MB-classifier gate only) -------------------
            if (needZdcRawForMinBias)
            {
                auto* inZdc = new Fun4AllNoSyncDstInputManager("DST_ZDC_RAW_IN");
                for (const auto& f : filesZdc) inZdc->AddFile(f);
                se->registerInputManager(inZdc);

                if (vlevel > 0)
                {
                    std::cout << "[INFO] AuAu MinimumBiasClassifier gate enabled: registered paired DST_ZDC_RAW input stream"
                              << " (nFiles=" << filesZdc.size() << ")" << std::endl;
                }
            }
        }
    
    if (isSim)
    {
        // For SIM we REQUIRE a reco-vertex stream from DST_GLOBAL.
        // DST_MBD_EPD is OPTIONAL: if absent (or 'NONE'), downstream code
        // will fall back to GlobalVertexMap for the reco vertex.
        // G4Hits is OPTIONAL: if missing, we skip truth-photon matching QA.
        if (!listHasGlobal)
        {
            std::ostringstream os;
            os << "isSim requires a reco-vertex DST_GLOBAL stream paired 1:1 with calo files.\n"
            << "Input list must include these columns per line:\n"
            << "  <DST_CALO_CLUSTER> <G4_or_NONE> <DST_JETS_or_NONE> <DST_GLOBAL> [<DST_MBD_EPD_or_NONE>]\n"
            << "But your list is missing DST_GLOBAL on at least one line.";
            detail::bail(os.str());
        }
        
        // G4 is OPTIONAL (photonjet productions may not provide it).
        // If you want to force it: export RJ_REQUIRE_G4=1
        bool requireG4 = false;
        if (const char* env = std::getenv("RJ_REQUIRE_G4")) requireG4 = (std::atoi(env) != 0);
        
        if (listHasG4 && (!auauCandidateSkimOnly || candidateSkimNeedsPhotonStitchTruth))
        {
            auto* inG4 = isSimEmbedded
            ? static_cast<Fun4AllInputManager*>(new Fun4AllNoSyncDstInputManager("DST_G4HITS_IN"))
            : static_cast<Fun4AllInputManager*>(new Fun4AllDstInputManager("DST_G4HITS_IN"));
            for (const auto& f : filesG4) inG4->AddFile(f);
            se->registerInputManager(inG4);
        }
        else if (listHasG4 && auauCandidateSkimOnly)
        {
            std::cout << "[INFO] RJ_AUAU_CANDIDATE_SKIM_ONLY=1: "
                      << "skipping the G4Hits input stream for the embedded-inclusive "
                      << "lane; truth-jet ownership is supplied by DST_JETS.\n";
        }
        else
        {
            if (requireG4 && !auauCandidateSkimOnly)
            {
                detail::bail("RJ_REQUIRE_G4=1 but no G4Hits stream was provided in the input list.");
            }
            else
            {
                std::cout << "\033[31m[WARN] isSim: no G4Hits stream provided (or it's 'NONE'). "
                "Continuing; truth-photon matching QA will be skipped.\033[0m\n";
            }
        }
        
        auto* inGlobal = isSimEmbedded
        ? static_cast<Fun4AllInputManager*>(new Fun4AllNoSyncDstInputManager("DST_GLOBAL_IN"))
        : static_cast<Fun4AllInputManager*>(new Fun4AllDstInputManager("DST_GLOBAL_IN"));
        for (const auto& f : filesGlobal) inGlobal->AddFile(f);
        se->registerInputManager(inGlobal);
        
        if (listHasMbd)
        {
            auto* inMbd = isSimEmbedded
            ? static_cast<Fun4AllInputManager*>(new Fun4AllNoSyncDstInputManager("DST_MBD_EPD_IN"))
            : static_cast<Fun4AllInputManager*>(new Fun4AllDstInputManager("DST_MBD_EPD_IN"));
            for (const auto& f : filesMbd) inMbd->AddFile(f);
            se->registerInputManager(inMbd);
        }
        else
        {
            std::cout << "\033[31m[WARN] isSim: no DST_MBD_EPD stream provided (or it's 'NONE'). "
            "Continuing; reco vertex will use GlobalVertexMap fallback.\033[0m\n";
        }
        
        if (verbose)
            std::cout << "[INFO] isSim: registered input managers (Calo + Global"
            << (listHasMbd ? " + MBD_EPD" : " (no MBD_EPD)")
            << ((listHasG4 && (!auauCandidateSkimOnly || candidateSkimNeedsPhotonStitchTruth))
                    ? " + G4" : " (no G4)") << ")\n";
    }
    
    
    // ------------------ Jets DST (SIM optional; for truth jets) -------
        if (isSim && useDSTTruthJets)
        {
            if (!listHasJets)
            {
                std::ostringstream os;
            os << "RJ_TRUTH_JETS_MODE=" << truthMode << " requires a list with at least 3 columns:\n"
            << "  <DST_CALO_CLUSTER> <G4Hits> <DST_JETS> [<DST_GLOBAL> <DST_MBD_EPD>]\n"
            << "Use your staged 5-column master list built from the matched lists.";
            detail::bail(os.str());
        }
        
        auto* inJets = isSimEmbedded
        ? static_cast<Fun4AllInputManager*>(new Fun4AllNoSyncDstInputManager("DST_JETS_IN"))
        : static_cast<Fun4AllInputManager*>(new Fun4AllDstInputManager("DST_JETS_IN"));
        for (const auto& f : filesJets) inJets->AddFile(f);
        se->registerInputManager(inJets);
        
        if (verbose)
            std::cout << "[INFO] isSim: registered DST_JETS input manager (truth jets from DST)\n";
    }
        }
    }
    
    if (verbose && isSim)
    {
        std::cout << "[INFO] Truth-jet mode (RJ_TRUTH_JETS_MODE=" << truthMode << "): "
        << (useDSTTruthJets ? "DST" : "")
        << ((useDSTTruthJets && buildTruthJetsFromParticles) ? "+" : "")
        << (buildTruthJetsFromParticles ? "BUILD" : "")
        << "\n";
    }

#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
    const std::string scaledTriggerOnlyEnv = env_lower("RJ_SCALED_TRIGGER_STUDY_ONLY");
    const std::string scaledTriggerCentEnv = env_lower("RJ_SCALED_TRIGGER_CENT_STUDY");
    const bool scaledTriggerCentStudy =
        (scaledTriggerCentEnv == "1" || scaledTriggerCentEnv == "true" ||
         scaledTriggerCentEnv == "yes" || scaledTriggerCentEnv == "on");
    const bool scaledTriggerStudyOnly =
        (scaledTriggerOnlyEnv == "1" || scaledTriggerOnlyEnv == "true" ||
         scaledTriggerOnlyEnv == "yes" || scaledTriggerOnlyEnv == "on" ||
         scaledTriggerCentStudy);
    if (scaledTriggerStudyOnly && (!isAuAuRequested || isSim))
    {
        detail::bail("RJ_SCALED_TRIGGER_STUDY_ONLY/RJ_SCALED_TRIGGER_CENT_STUDY is valid only for AuAu data.");
    }
#endif
    
    // --------------------------------------------------------------------
    // 3.  Geometry + status + calibration + clustering
    // --------------------------------------------------------------------
    //
    // IMPORTANT:
    // RawClusterBuilderTemplate::process_event() ALWAYS requires "TOWERGEOM_CEMC".
    // But CaloGeomMapping with UseDetailedGeometry(true) publishes ONLY
    // "TOWERGEOM_CEMC_DETAILED" for CEMC.
    // So we ALWAYS create the legacy node (TOWERGEOM_CEMC), and optionally
    // also create the detailed node (TOWERGEOM_CEMC_DETAILED).
    //
    if (!ppg12ClosureCanary)
    {
        bool useDetailedCemcGeom = true;  // keep ordinary production behavior
        if (const char* env = std::getenv("RJ_DETAILED_CEMC_GEOM"))
        {
            useDetailedCemcGeom = (std::atoi(env) != 0);
        }

    // The deployed PPG12 run-28 DI chain lets Process_Calo_Calib publish its
    // own geometry.  Registering the unified helpers first changes that graph.
    if (!usePPG12ArchivedDIG4OnlyReco)
    {
        // Always publish the legacy/simple CEMC geometry node: TOWERGEOM_CEMC
        {
            auto* geomCemcLegacy = new CaloGeomMapping("Geom_CEMC");
            geomCemcLegacy->set_detector_name("CEMC");
            geomCemcLegacy->set_UseDetailedGeometry(false);
            se->registerSubsystem(geomCemcLegacy);
        }

        // Optionally publish the detailed CEMC geometry node.
        if (useDetailedCemcGeom)
        {
            auto* geomCemcDetailed = new CaloGeomMapping("Geom_CEMC_DETAILED");
            geomCemcDetailed->set_detector_name("CEMC");
            geomCemcDetailed->set_UseDetailedGeometry(true);
            se->registerSubsystem(geomCemcDetailed);
        }

        // HCAL nodes (detailed not supported; mapping falls back internally).
        for (const std::string& det : {"HCALIN","HCALOUT"})
        {
            auto* geom = new CaloGeomMapping(("Geom_" + det).c_str());
            geom->set_detector_name(det);
            geom->set_UseDetailedGeometry(true);
            se->registerSubsystem(geom);
        }
    }
    else if (verbose || vlevel > 0)
    {
        std::cout << "[PPG12_DI_ARCHIVED_RECO] Process_Calo_Calib owns the "
                  << "deployed run-28 geometry contract\n";
        }
    }
    
    
    // ------------------------------------------------------------------
    // Calo calibration + clustering contract
    //
    //   jetcalo     -> lower-level tower input: run Process_Calo_Calib()
    //                  (also the default for embedded DST_CALO MC)
    //   calofitting -> waveform-fit input: run Process_Calo_Calib()
    //   simdst      -> analysis DST already carries calibrated towers/clusters
    //   RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
    //              -> pp SIM diagnostic path: reproduce PPG12 anatreemaker's
    //                 G4/input/pedestal stack and rebuild CEMC clusters.
    // ------------------------------------------------------------------
    if (usePPG12PPSimRebuildCaloFromG4)
    {
        if (vlevel > 0)
        {
            std::cout << "[PPG12_PPSIM_REBUILD_CALO_FROM_G4] "
                      << ((ppg12ClosureCanary || !usePPG12Fig11G4OnlyRebuild)
                              ? (usePPG12ArchivedDIG4OnlyReco
                                     ? "running Process_Calo_Calib() on the deployed PPG12 DI graph\n"
                                     : "running Process_Calo_Calib() after the PPG12 input/pedestal stack\n")
                              : "skipping Process_Calo_Calib(); the non-canary Fig.11 G4-only path already registered its tower and cluster stack\n");
        }
#if RJ_HAS_SPHENIX_G4_INPUT_MACROS
        if (ppg12ClosureCanary || !usePPG12Fig11G4OnlyRebuild)
        {
            Process_Calo_Calib();
        }

        // Only the exact executable-oracle branch adds the preserved
        // post-calibration no-split collection.  Ordinary G4-rebuild jobs
        // retain their pre-canary reconstruction graph.
        if (ppg12ClosureCanary)
        {
            auto* noSplitBuilder =
                new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate_PPG12OracleNoSplit");
            noSplitBuilder->Detector("CEMC");
            noSplitBuilder->set_threshold_energy(0.070);
            const char* calibrationRoot = std::getenv("CALIBRATIONROOT");
            if (!calibrationRoot || !std::string(calibrationRoot).size())
                detail::bail("CALIBRATIONROOT is required for the PPG12 no-split cluster builder");
            const std::string emcProfile =
                std::string(calibrationRoot) +
                "/EmcProfile/CEMCprof_Thresh30MeV.root";
            noSplitBuilder->LoadProfile(emcProfile);
            noSplitBuilder->setSubclusterSplitting(false);
            noSplitBuilder->setOutputClusterNodeName("CLUSTERINFO_CEMC_NO_SPLIT");
            noSplitBuilder->set_UseTowerInfo(1);
            noSplitBuilder->Verbosity(0);
            se->registerSubsystem(noSplitBuilder);

            if (verbose || vlevel > 0)
            {
                std::cout << "[PPG12_PPSIM_REBUILD_CALO_FROM_G4] registered "
                          << "post-calibration no-split TowerInfo cluster builder "
                          << "node=CLUSTERINFO_CEMC_NO_SPLIT" << std::endl;
            }
        }
#else
        detail::bail(
            "RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1 requested tower rebuild, but "
            "GlobalVariables.C/G4_Input.C were not available to this ROOT macro environment.");
#endif
    }
    else if (isSim && caloInputMode == "simdst")
    {
        if (vlevel > 0)
        {
            std::cout << (isSimEmbedded ? "[isSimEmbedded]" : "[isSim]")
            << " skipping Process_Calo_Calib() "
            << "(SIM DST already has calibrated towers/clusters)\n";
        }
    }
    else if (caloInputMode == "calofitting" || caloInputMode == "jetcalo")
    {
        if (vlevel > 0)
        {
            if (isSimEmbedded)
            {
                if (caloInputMode == "calofitting")
                {
                    std::cout << "[isSimEmbedded][embedded DST_CALO contract] "
                              << "running Process_Calo_Calib() on CALOFITTING input\n";
                }
                else
                {
                    std::cout << "[isSimEmbedded][embedded DST_CALO contract] "
                              << "running Process_Calo_Calib() on DST_CALO / JETCALO input\n";
                }
            }
            else if (caloInputMode == "calofitting")
            {
                std::cout << "[DATA] running Process_Calo_Calib() on CALOFITTING input\n";
            }
            else
            {
                std::cout << "[DATA] running Process_Calo_Calib() on JETCALO input\n";
            }
            if (isAuAuRequested)
            {
                std::cout << "[DATA][AuAu] clusterUEpipeline may still apply native UE subtraction and reclusterization afterward\n";
            }
        }
        if (!isSim && isAuAuRequested && caloInputMode == "calofitting")
        {
            setenv("RJ_SKIP_CALO_STATUS_SKIMMER", "1", 1);
            unsetenv("RJ_SKIP_CALO_TOWER_STATUS");
            unsetenv("RJ_DEFER_CALO_TOWER_STATUS_TO_CALIB");
            unsetenv("RJ_USE_TOWERINFO_CALO_INPUT");
            unsetenv("RJ_CALO_INPUT_NODE_PREFIX");
            setenv("RJ_DISABLE_CEMC_BAD_TOWER_MASK", "1", 1);
            if (vlevel > 0)
            {
                std::cout << "[DATA][AuAu] CALOFITTING/TowerInfo input: skipping legacy CaloStatusSkimmer, "
                          << "running the stock CaloTowerStatus -> CaloTowerCalib chain on TOWERS_* so status is copied "
                          << "into TOWERINFO_CALIB_* before cluster building, and disabling the downstream "
                          << "PhotonClusterBuilder CEMC mask for canonical proof\n";
            }
        }
        if (!isSim && (usePPG12PPDataPair || useAuAuJetCaloDataPair) && caloInputMode == "jetcalo")
        {
            // The paired DST_JET + DST_JETCALO data streams already carry
            // TOWERINFO_CALIB_* nodes. Keep those producer status/calibration
            // nodes intact instead of rerunning the legacy status setter on
            // missing raw TOWERS_* input.
            setenv("RJ_SKIP_CALO_TOWER_STATUS", "1", 1);
            unsetenv("RJ_CALO_TOWER_STATUS_INPUT_PREFIX");
            if (vlevel > 0)
            {
                std::cout << (useAuAuJetCaloDataPair ? "[DATA][AuAu]" : "[DATA][PPG12 pp]")
                          << " JETCALO/TowerInfo input: "
                          << "preserving existing TOWERINFO_CALIB_* status/calibration nodes "
                          << (useAuAuJetCaloDataPair ? "(DST_JET/DST_JETCALO data policy)\n" : "(PPG12 ana521 parity)\n");
            }
        }
        if (isSimEmbedded)
        {
            setenv("RJ_DISABLE_CEMC_BAD_TOWER_MASK", "1", 1);
            if (vlevel > 0)
            {
                std::cout << "[isSimEmbedded] disabling downstream PhotonClusterBuilder CEMC bad-tower mask; "
                          << "embedded MC uses the producer tower-quality state while Process_Calo_Calib "
                          << "provides the calibrated CEMC/HCAL nodes needed by clustering and RetowerCEMC\n";
            }
        }
        std::shared_ptr<rj_cemc_status::State> cemcStatusState;
        if (useAuAuJetCaloDataPair)
        {
            cemcStatusState = std::make_shared<rj_cemc_status::State>();
            cemcStatusState->expected_run = run;
            if (filesZdc.size() != 1)
                detail::bail("AuAu status-validated capture requires one exact DST_JET/DST_JETCALO pair per invocation");
            cemcStatusState->original_calo_input = filesZdc.front();
            cemcStatusState->recovery_authorized = env_truthy_local("RJ_AUAU_CEMC_RESTORE_MISSING_STATUS");
            const char* recoveryMap = std::getenv("RJ_AUAU_CEMC_RECOVERY_MAP");
            cemcStatusState->recovery_payload = recoveryMap ? recoveryMap : "";
            if (cemcStatusState->recovery_authorized != !cemcStatusState->recovery_payload.empty())
                detail::bail("AuAu CEMC recovery requires both explicit repair mode and a pinned run-specific payload; no implicit current-tag fallback");
        }
        Process_Calo_Calib(cemcStatusState);
    }
    else
    {
        detail::bail("unsupported caloInputMode '" + caloInputMode + "'");
    }
    
    
    if (isSimEmbedded)
    {
        if (vlevel > 0)
        {
            std::cout << "[isSimEmbedded] skipping MbdReco (use embedded sample's existing MBD products)" << std::endl;
            std::cout << "[isSimEmbedded] skipping ZdcReco (not needed for embedded minimal path)" << std::endl;
            std::cout << "[isSimEmbedded] skipping GlobalVertexReco (use embedded sample's existing GlobalVertexMap)" << std::endl;
        }
    }
    else if (usePPG12PPSimRebuildCaloFromG4)
    {
        if (vlevel > 0)
        {
            std::cout << "[PPG12_PPSIM_REBUILD_CALO_FROM_G4] MBD/GlobalVertex reconstruction already registered early "
                      << (usePPG12PPSimG4OnlyInput ? "for G4-only pp SIM input" : "for four-lane PPG12 pp SIM input")
                      << std::endl;
        }
    }
    else
    {
        if (vlevel > 0) std::cout << "Calibrating MBD" << std::endl;
        std::unique_ptr<MbdReco> mbdreco = std::make_unique<MbdReco>();
        se->registerSubsystem(mbdreco.release());
        
        if (!isSim && !isPPrun25)
        {
            if (vlevel > 0)
            {
                if (needZdcRawForMinBias)
                    std::cout << "Calibrating ZDC after Process_Calo_Calib ZDC tower building" << std::endl;
                else
                    std::cout << "Calibrating ZDC" << std::endl;
            }

            auto* zdcreco = new ZdcReco();
            zdcreco->set_zdc1_cut(0.0);
            zdcreco->set_zdc2_cut(0.0);
            se->registerSubsystem(zdcreco);
        }
        else
        {
            if (vlevel > 0)
            {
                if (isPPrun25) std::cout << "[isPPrun25] skipping ZdcReco (CALOFITTING DST may not have TOWERS_ZDC)" << std::endl;
                else           std::cout << "[isSim] skipping ZdcReco (sim DST has no TOWERS_ZDC)" << std::endl;
            }
        }
        
        
        if (vlevel > 0) std::cout << "Retrieving Vtx Info" << std::endl;
        std::unique_ptr<GlobalVertexReco> gvertex = std::make_unique<GlobalVertexReco>();
        se->registerSubsystem(gvertex.release());
    }
    
    bool isAuAuData = false;
    bool isAuAuLike = false;
    if (isSimEmbedded)
    {
        isAuAuData = false;
        isAuAuLike = true;
    }
    else if (isSim)
    {
        isAuAuData = false;
        isAuAuLike = false;
    }
    else if (const char* env = std::getenv("RJ_DATASET"))
    {
        std::string sLower = detail::trim(std::string(env));
        std::transform(sLower.begin(), sLower.end(), sLower.begin(), [](unsigned char c){ return std::tolower(c); });
        if (sLower == "isauau" || sLower == "auau" || sLower == "aa")
            isAuAuData = true;
        else if (sLower == "ispp" || sLower == "pp" || sLower == "ispprun25" || sLower == "pprun25" || sLower == "pp25")
            isAuAuData = false;
        else
            isAuAuData = (run > 53864);
        isAuAuLike = isAuAuData;
    }
    else
    {
        isAuAuData = (run > 53864);
        isAuAuLike = isAuAuData;
    }
    
    if (verbose || vlevel > 0)
    {
        std::cout << "[FLOW] dataset semantics:"
        << " | isSim=" << (isSim ? "true" : "false")
        << " | isSimEmbedded=" << (isSimEmbedded ? "true" : "false")
        << " | isAuAuData=" << (isAuAuData ? "true" : "false")
        << " | isAuAuLike=" << (isAuAuLike ? "true" : "false")
        << " | isPPrun25=" << (isPPrun25 ? "true" : "false")
        << std::endl;
    }

#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
    const bool centralityWillRun =
        isAuAuData && (!scaledTriggerStudyOnly || scaledTriggerCentStudy);
#else
    const bool centralityWillRun = false;
#endif
    const auto centralityBinding =
        rj_centrality_binding::resolve(isAuAuData, centralityWillRun, run);

#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
    if (scaledTriggerStudyOnly)
    {
        if (vlevel > 0)
        {
            std::cout << "[FLOW] scaled-trigger study only: registering lightweight max-cluster module\n"
                      << "       keeps geometry, calo calibration/clustering, MBD reco, and global vertex reco\n"
                      << "       " << (scaledTriggerCentStudy ? "builds MinimumBiasClassifier/CentralityReco for centrality-sliced trigger QA"
                                                           : "skips centrality")
                      << ", photon selection, UE subtraction, jet reco, recoil analysis, and trees\n";
        }

        if (scaledTriggerCentStudy)
        {
            if (vlevel > 0)
            {
                std::cout << "[FLOW] scaled-trigger centrality study: registering MinimumBiasClassifier and CentralityReco\n";
            }
            auto* mb = new MinimumBiasClassifier();
            mb->Verbosity(0);
            se->registerSubsystem(mb);

            auto* cent = new CentralityReco();
            cent->Verbosity(0);
            rj_centrality_binding::apply(cent, centralityBinding);
            se->registerSubsystem(cent);
        }

        se->registerSubsystem(new ScaledTriggerStudyReco(outRoot));

        try
        {
            if (vlevel > 0) std::cout << "[INFO] Starting scaled-trigger-only event loop ..." << std::endl;
            const int runRc = se->run(nEvents);
            // Keep the event loop quiet in batch, but never suppress End()
            // diagnostics. Fun4AllServer::End() returns only an aggregate, so
            // subsystem output is required to classify a nonzero value safely.
            _silence.disable();
            if (vlevel > 0) std::cout << "[INFO] Calling se->End() ..." << std::endl;
            const int endRc = se->End();
            detail::enforce_fun4all_status(
                "scaled-trigger-only",
                se,
                nEvents,
                runRc,
                endRc,
                permittedRepeatingPedestalInputManager);
            if (vlevel > 0) std::cout << "[INFO] Finished scaled-trigger-only job." << std::endl;
        }
        catch (const std::exception& e)
        {
            detail::bail(std::string("exception in scaled-trigger-only Fun4All path: ") + e.what());
        }
        return;
    }
#endif
    
    if (isAuAuLike && !isSimEmbedded)
    {
        if (cfg.setMinBiasClassifer)
        {
            if (vlevel > 0) std::cout << "building minbias classifier" << std::endl;
            auto* mb = new MinimumBiasClassifier();
            mb->Verbosity(0);
            se->registerSubsystem(mb);
        }
        else if (vlevel > 0)
        {
            std::cout << "[AuAu] skipping MinimumBiasClassifier gate (setMinBiasClassifer=false)" << std::endl;
        }
        
        if (vlevel > 0) std::cout << "building centrality classifier (Au+Au-like)" << std::endl;
        auto* cent = new CentralityReco();
        cent->Verbosity(0);
        rj_centrality_binding::apply(cent, centralityBinding);
        se->registerSubsystem(cent);
    }
    else
    {
        if (vlevel > 0)
        {
            if (isSimEmbedded) std::cout << "[isSimEmbedded] skipping MinimumBiasClassifier/CentralityReco (use embedded sample's existing centrality products)" << std::endl;
            else               std::cout << "[pp dataset] skipping CentralityReco" << std::endl;
        }
    }
    
    setenv("BEMCREC_CEMC_DISABLE_ASINH_POSITION", "0", 1);
    
    if (cfg.doPi0Analysis)
    {
        if (vlevel > 0)
        {
            std::cout << "[pi0] position-corrected-only mode: using CLUSTERINFO_CEMC only"
                      << " (no parallel CLUSTERINFO_CEMC_NOCORR branch will be built)" << std::endl;
        }
    }
    
    // ---------------------- Reco jets -----------------------------------------
    // Au+Au:
    //   - tower-level UE subtraction (RetowerCEMC + DTB + CASJ + DTB2 + SubtractTowers)
    //   - final jets are built from SUB1 tower containers
    //   - jets are written to: AntiKt_Tower_<rKey>_Sub1  (e.g. AntiKt_Tower_r04_Sub1)
    //
    // PP / isSim:
    //   - keep the existing pp-style JetReco + JetCalib chain (see 'else' below)
    // --------------------------------------------------------------------------
    std::string towerPrefixPCB = "TOWERINFO_CALIB";
    if (isAuAuLike)
    {
        if (const char* env = std::getenv("RJ_TOWERINFO_PREFIX"))
        {
            std::string s = detail::trim(std::string(env));
            if (!s.empty()) towerPrefixPCB = s;
        }
    }
    
    if (verbose || vlevel > 0)
    {
        std::cout << "[FLOW] reco/calibration branch:"
        << " | branch=" << (isAuAuLike ? "AuAu-like HI UE subtraction + SUB1 jets + JetCalib"
                            : "pp-style jets + JetCalib")
        << " | Process_Calo_Calib=" << (usePPG12PPSimRebuildCaloFromG4 ||
                                       (!isSim && (caloInputMode == "calofitting" || caloInputMode == "jetcalo")) ||
                                       (isSimEmbedded && (caloInputMode == "calofitting" || caloInputMode == "jetcalo"))
                                           ? "ON" : "OFF")
        << " | clusterUEpipeline=" << cfg.clusterUEpipeline
        << " | towerPrefixPCB=" << towerPrefixPCB
        << " | truthJets=" << (useDSTTruthJets ? "DST" : "BUILD")
        << ((useDSTTruthJets && buildTruthJetsFromParticles) ? "+BUILD" : "")
        << std::endl;
    }
    
    if (isAuAuLike)
    {
        // Ensure JetBackground modules can attach nodes under DST/TOWER
        auto* ensure = new EnsureJetCalibNodes("EnsureJetCalibNodes_forHIJets");
        ensure->Verbosity(0);
        se->registerSubsystem(ensure);
        
        // Optional: control HI UE-subtraction verbosity independently
        int hiV = 0;
        if (const char* env = std::getenv("RJ_HIUE_VERBOSITY")) hiV = std::atoi(env);
        
        // Match Macro_HIJetReco.C semantics:
        //   0 = no flow
        //   1 = psi2 derived from calo
        //   2 = psi2 derived from HIJING
        //   3 = psi2 derived from sEPD
        int hiFlow = 0;
        if (const char* env = std::getenv("RJ_HI_DO_FLOW")) hiFlow = std::atoi(env);
        
        // TowerInfo node prefix for HI background chain.
        // calo-fitting Au+Au DSTs often publish per-detector nodes as:
        //   TOWERINFO_CEMC, TOWERINFO_HCALIN, TOWERINFO_HCALOUT
        // (NOT TOWERINFO_CALIB_*)
        std::string towerPrefix = "TOWERINFO_CALIB";
        if (const char* env = std::getenv("RJ_TOWERINFO_PREFIX"))
        {
            std::string s = detail::trim(std::string(env));
            if (!s.empty()) towerPrefix = s;
        }
        
        if (vlevel > 0)
            std::cout << "[HI] UE subtraction enabled: towerPrefix=" << towerPrefix
            << " do_flow=" << hiFlow
            << " (HIUE Verbosity=" << hiV << ")\n";
        
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
        if (hiFlow == 3)
        {
            auto* epreco = new EventPlaneReco();
            se->registerSubsystem(epreco);
        }
#endif
        
        // ------------------------------------------------------------------
        // 1) Retower CEMC (towerinfo)
        // ------------------------------------------------------------------
        auto* rcemc = new RetowerCEMC();
        rcemc->Verbosity(hiV);
        rcemc->set_towerinfo(true);
        rcemc->set_frac_cut(0.5);
        rcemc->set_towerNodePrefix(towerPrefix);
        se->registerSubsystem(rcemc);
        
        // ------------------------------------------------------------------
        // 2) Seed jets (RAW, R=0.2) for background estimation
        //     -> AntiKt_TowerInfo_HIRecoSeedsRaw_r02
        // ------------------------------------------------------------------
        {
            auto* seedReco = new JetReco("JetsReco_HIRecoSeedsRaw_r02");
            
            auto* incemc  = new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER, towerPrefix);
            auto* inihcal = new TowerJetInput(Jet::HCALIN_TOWERINFO,       towerPrefix);
            auto* inohcal = new TowerJetInput(Jet::HCALOUT_TOWERINFO,      towerPrefix);
            
            incemc->set_GlobalVertexType(GlobalVertex::MBD);
            inihcal->set_GlobalVertexType(GlobalVertex::MBD);
            inohcal->set_GlobalVertexType(GlobalVertex::MBD);
            
            seedReco->add_input(incemc);
            seedReco->add_input(inihcal);
            seedReco->add_input(inohcal);
            
            seedReco->add_algo(detail::fjAlgo(0.2f), "AntiKt_TowerInfo_HIRecoSeedsRaw_r02");
            seedReco->set_algo_node("ANTIKT");
            seedReco->set_input_node("TOWER");
            seedReco->Verbosity(hiV);
            se->registerSubsystem(seedReco);
        }
        
        // ------------------------------------------------------------------
        // 3) Background iteration 1 (seedType=0 -> RAW seeds)
        //     -> TowerInfoBackground_Sub1
        // ------------------------------------------------------------------
        auto* dtb = new DetermineTowerBackground();
        dtb->SetBackgroundOutputName("TowerInfoBackground_Sub1");
        dtb->SetFlow(hiFlow);
        dtb->SetSeedType(0);
        dtb->SetSeedJetD(3);
        dtb->Verbosity(hiV);
        dtb->set_towerNodePrefix(towerPrefix);
        se->registerSubsystem(dtb);
        
        // ------------------------------------------------------------------
        // 4) Copy + subtract jets using Sub1 background
        //     -> creates AntiKt_TowerInfo_HIRecoSeedsSub_r02
        // ------------------------------------------------------------------
        auto* casj = new CopyAndSubtractJets();
        casj->SetFlowModulation(hiFlow);
        casj->Verbosity(hiV);
        casj->set_towerinfo(true);
        casj->set_towerNodePrefix(towerPrefix);
        se->registerSubsystem(casj);
        
        // ------------------------------------------------------------------
        // 5) Background iteration 2 (seedType=1 -> SUB seeds)
        //     -> TowerInfoBackground_Sub2
        // ------------------------------------------------------------------
        auto* dtb2 = new DetermineTowerBackground();
        dtb2->SetBackgroundOutputName("TowerInfoBackground_Sub2");
        dtb2->SetFlow(hiFlow);
        dtb2->SetSeedType(1);
        dtb2->SetSeedJetPt(7);
        dtb2->Verbosity(hiV);
        dtb2->set_towerNodePrefix(towerPrefix);
        se->registerSubsystem(dtb2);
        
        // ------------------------------------------------------------------
        // 6) Subtract towers using Sub2 background
        //     -> writes *SUB1 tower containers used for final jet reco
        // ------------------------------------------------------------------
        auto* st = new SubtractTowers();
        st->SetFlowModulation(hiFlow);
        st->Verbosity(hiV);
        st->set_towerinfo(true);
        st->set_towerNodePrefix(towerPrefix);
        se->registerSubsystem(st);
        
        // ------------------------------------------------------------------
        // 7) Parallel final jet views (one node pair per R and view)
        //
        //     Un-subtracted retowered towers:
        //       AntiKt_Tower_<rKey>_NoSub_RAW  (before CDB JES)
        //       AntiKt_Tower_<rKey>_NoSub      (after CDB JES)
        //
        //     SUB1 UE-subtracted towers:
        //     -> RAW jets written to AntiKt_Tower_<rKey>_Sub1_RAW
        //     -> JetCalib output written to AntiKt_Tower_<rKey>_Sub1
        // ------------------------------------------------------------------
        int jetcalV = 0;
        if (const char* env = std::getenv("RJ_JETCALIB_VERBOSITY")) jetcalV = std::atoi(env);

        if (auauCandidateSkimOnly && vlevel > 0)
        {
            std::cout << "[THE-105 candidate skim] retaining AuAu tower subtraction "
                      << "for the photon builder; final reconstructed jets and "
                      << "JetCalib are disabled" << std::endl;
        }

        if (!auauCandidateSkimOnly)
        for (const auto& radKey : activeJetRKeys)
        {
            int D = 0;
            try { D = std::stoi(radKey.substr(1)); } catch (...) { continue; }
            if (D <= 0) continue;
            const float R = 0.1f * D;

            const std::string calibNode =
                std::string("AntiKt_Tower_") + radKey + "_NoSub" +
                analysisJetNodeSuffix;
            const std::string rawNode = calibNode + "_RAW";
            const std::string recoName =
                std::string("JetsReco_AuAuNoSub_") + radKey;

            auto* jreco = new JetReco(recoName);
            auto* incemc =
                new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER,towerPrefix);
            auto* inihcal =
                new TowerJetInput(Jet::HCALIN_TOWERINFO,towerPrefix);
            auto* inohcal =
                new TowerJetInput(Jet::HCALOUT_TOWERINFO,towerPrefix);
            incemc->set_GlobalVertexType(GlobalVertex::MBD);
            inihcal->set_GlobalVertexType(GlobalVertex::MBD);
            inohcal->set_GlobalVertexType(GlobalVertex::MBD);
            jreco->add_input(incemc);
            jreco->add_input(inihcal);
            jreco->add_input(inohcal);
            jreco->add_algo(detail::fjAlgo(R, true),rawNode);
            jreco->set_algo_node("ANTIKT");
            jreco->set_input_node("TOWER");
            int jetrecoV = 0;
            if (const char* env = std::getenv("RJ_JETRECO_VERBOSITY"))
                jetrecoV = std::atoi(env);
            jreco->Verbosity(jetrecoV);
            se->registerSubsystem(jreco);

            auto* jcal =
                new JetCalib(std::string("JetCalib_AuAuNoSub_") + radKey);
            jcal->set_InputNode(rawNode);
            jcal->set_OutputNode(calibNode);
            jcal->set_JetRadius(R);
            jcal->set_ApplyZvrtxDependentCalib(true);
            jcal->set_ApplyEtaDependentCalib(true);
            jcal->Verbosity(jetcalV);
            se->registerSubsystem(jcal);

            auto* probe = new JetCalibOneEventProbe(
                std::string("JetCalibOneEventProbe_AuAuNoSub_") + radKey,
                rawNode,calibNode,/*maxJetsToPrint=*/12);
            probe->Verbosity(vlevel);
            se->registerSubsystem(probe);

            if (vlevel > 0)
                std::cout << "[INFO] (AuAu) reco jets: built " << rawNode
                          << " -> " << calibNode << " (R=" << R
                          << ") from un-subtracted retowered towers with JetCalib\n";
        }

        if (!auauCandidateSkimOnly)
        for (const auto& radKey : activeJetRKeys)
        {
            int D = 0;
            try { D = std::stoi(radKey.substr(1)); } catch (...) { continue; }
            if (D <= 0) continue;
            const float R = 0.1f * D;
            
            const std::string calibNode = std::string("AntiKt_Tower_") + radKey + "_Sub1" + analysisJetNodeSuffix;
            const std::string rawNode   = calibNode + "_RAW";
            const std::string recoName  = std::string("JetsReco_AuAuSub_") + radKey;
            
            auto* jreco = new JetReco(recoName);
            
            auto* incemc  = new TowerJetInput(Jet::CEMC_TOWERINFO_SUB1,   towerPrefix);
            auto* inihcal = new TowerJetInput(Jet::HCALIN_TOWERINFO_SUB1, towerPrefix);
            auto* inohcal = new TowerJetInput(Jet::HCALOUT_TOWERINFO_SUB1, towerPrefix);
            
            incemc->set_GlobalVertexType(GlobalVertex::MBD);
            inihcal->set_GlobalVertexType(GlobalVertex::MBD);
            inohcal->set_GlobalVertexType(GlobalVertex::MBD);
            
            jreco->add_input(incemc);
            jreco->add_input(inihcal);
            jreco->add_input(inohcal);
            
            jreco->add_algo(detail::fjAlgo(R, true), rawNode);
            jreco->set_algo_node("ANTIKT");
            jreco->set_input_node("TOWER");
            
            int jetrecoV = 0;
            if (const char* env = std::getenv("RJ_JETRECO_VERBOSITY")) jetrecoV = std::atoi(env);
            jreco->Verbosity(jetrecoV);
            
            se->registerSubsystem(jreco);
            
            {
                auto* jcal = new JetCalib(std::string("JetCalib_AuAuSub_") + radKey);
                jcal->set_InputNode(rawNode);
                jcal->set_OutputNode(calibNode);
                jcal->set_JetRadius(R);
                jcal->set_ApplyZvrtxDependentCalib(true);
                jcal->set_ApplyEtaDependentCalib(true);
                jcal->Verbosity(jetcalV);
                se->registerSubsystem(jcal);
                
                auto* probe = new JetCalibOneEventProbe(std::string("JetCalibOneEventProbe_AuAuSub_") + radKey,
                                                        rawNode,
                                                        calibNode,
                                                        /*maxJetsToPrint=*/12);
                probe->Verbosity(vlevel);
                se->registerSubsystem(probe);
                
                if (vlevel > 0)
                    std::cout << "[INFO] (AuAu) reco jets: built " << rawNode << " -> " << calibNode << " (R=" << R << ") from SUB1 towers with JetCalib\n";
            }
        }
    }
    else
    {
        // ---------------------- Reco jets + JES calibration (pp-style only) ----------------------
        //
        // IMPORTANT:
        // JetCalib::CreateNodeTree() requires a PHCompositeNode named "TOWER".
        // Many pp DSTs do NOT have it, so JetCalib aborts unless we create it.
        // We register EnsureJetCalibNodes once (only when we intend to run JetCalib).
        //
        // Apply pp JES calibration for ALL pp-style running (pp data AND isSim).
        // Keep Au+Au-like chains excluded.
        const bool doJetCalibAny = (!isAuAuLike);
        if (doJetCalibAny)
        {
            auto* ensure = new EnsureJetCalibNodes("EnsureJetCalibNodes_forJES");
            // keep this modest; set RJ_JETCALIB_NODE_VERBOSE=1 for prints
            int nodeV = 0;
            if (const char* env = std::getenv("RJ_JETCALIB_NODE_VERBOSE")) nodeV = std::atoi(env);
            ensure->Verbosity(nodeV);
            se->registerSubsystem(ensure);
            
            if (vlevel > 0)
                std::cout << "[INFO] JES: enabling JetCalib (pp-style: pp data + isSim) => ensuring DST/TOWER exists\n";
        }
        
        // Optional: control JetCalib verbosity independently
        int jetcalV = 0; // default: show InitRun/process_event messages
        if (const char* env = std::getenv("RJ_JETCALIB_VERBOSITY"))
        {
            jetcalV = std::atoi(env);
        }
        
        for (const auto& radKey : activeJetRKeys)
        {
            int D = 0;
            try { D = std::stoi(radKey.substr(1)); } catch (...) { continue; }
            if (D <= 0) continue;
            const float R = 0.1f * D;
            
            // Canonical node name that RecoilJets reads (naming convention)
            const std::string calibNode = std::string("AntiKt_Tower_") + radKey + analysisJetNodeSuffix;
            
            // Apply JES calibration for pp-like chains (pp data + pp-style SIM)
            // Run pp JES calibration in BOTH pp data and isSim (pp-style chains)
            const bool doJetCalib = (!isAuAuLike);
            
            // If calibrating: build RAW jets to a separate node to avoid name collision
            const std::string rawNode = doJetCalib ? (calibNode + "_RAW") : calibNode;
            
            // ---------------------- JetReco (build jets) ----------------------
            const std::string recoName = std::string("JetsReco_") + (isSim ? "Sim_" : "Data_") + radKey;
            
            auto* jreco = new JetReco(recoName);
            jreco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO,    "TOWERINFO_CALIB"));
            jreco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO,  "TOWERINFO_CALIB"));
            jreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO, "TOWERINFO_CALIB"));
            
            jreco->add_algo(detail::fjAlgo(R, true), rawNode);
            jreco->set_algo_node("ANTIKT");
            jreco->set_input_node("TOWERINFO_CALIB");
            
            // If you want JetReco debug: export RJ_JETRECO_VERBOSITY=1
            int jetrecoV = 0;
            if (const char* env = std::getenv("RJ_JETRECO_VERBOSITY")) jetrecoV = std::atoi(env);
            jreco->Verbosity(jetrecoV);
            
            se->registerSubsystem(jreco);
            
            if (vlevel > 0)
            {
                std::cout << "[INFO] (" << (isSim ? "isSim" : "data") << ") reco jets: built "
                << rawNode << " (R=" << R << ") from TOWERINFO_CALIB\n";
            }
            
            // ---------------------- JetCalib (apply JES) ----------------------
            if (doJetCalib)
            {
                auto* jcal = new JetCalib(std::string("JetCalib_") + radKey);
                
                // JetCalib reads RAW jets and writes CALIB jets into the canonical node
                jcal->set_InputNode(rawNode);
                jcal->set_OutputNode(calibNode);
                
                jcal->set_JetRadius(R);
                
                // Full pp JES: Zvrtx + eta dependent
                jcal->set_ApplyZvrtxDependentCalib(true);
                jcal->set_ApplyEtaDependentCalib(true);
                
                jcal->Verbosity(jetcalV);
                se->registerSubsystem(jcal);
                
                // ------------------------------------------------------------
                // Targeted JES probe: prints ONCE, ONLY if Verbosity() >= 20
                // ------------------------------------------------------------
                auto* probe = new JetCalibOneEventProbe(std::string("JetCalibOneEventProbe_") + radKey,
                                                        rawNode,
                                                        calibNode,
                                                        /*maxJetsToPrint=*/12);
                probe->Verbosity(vlevel);  // uses your RJ_VERBOSITY; probe prints only if >=20
                se->registerSubsystem(probe);
                
                if (vlevel > 0)
                {
                    std::cout << "[INFO] JES calib enabled: " << rawNode << " -> " << calibNode
                    << " (R=" << R << ", Zvrtx+eta dependent, JetCalibVerbosity=" << jetcalV << ")\n";
                }
            }
        }
    }
    
    // ---------------------- Truth jets -----------------------------------------
    // If useDSTTruthJets==true: they already exist from DST_JETS_IN.
    // If buildTruthJetsFromParticles==true: build them from TRUTH particles here.
    if (isSim && buildTruthJetsFromParticles && !auauCandidateSkimOnly)
    {
        if (vlevel > 0)
        {
            if (useDSTTruthJets)
                std::cout << "[INFO] (isSim) truth jets: ALSO building from TRUTH particles (QA 'both' mode)\n";
            else
                std::cout << "[INFO] (isSim) truth jets: building from TRUTH particles (no DST_JETS)\n";
        }
        
        for (const auto& radKey : activeJetRKeys)
        {
            int D = 0;
            try { D = std::stoi(radKey.substr(1)); } catch (...) { continue; }
            if (D <= 0) continue;
            const float R = 0.1f * D;
            
            auto* truthReco = new JetReco(std::string("TruthJetReco_FromParticles_") + radKey);
            auto* tji = new TruthJetInput(Jet::PARTICLE);
            tji->add_embedding_flag(1);  // pythia/herwig particles only
            truthReco->add_input(tji);
            
            // Node naming:
            //  - If we're NOT reading DST jets, keep canonical name "AntiKt_Truth_<rKey>"
            //  - If we ARE reading DST jets too (BOTH mode), avoid collisions
            const std::string truthNode = (useDSTTruthJets && buildTruthJetsAsAltNode)
            ? (std::string("AntiKt_TruthFromParticles_") + radKey)
            : (std::string("AntiKt_Truth_") + radKey);
            
            truthReco->add_algo(detail::fjAlgo(R), truthNode);
            truthReco->set_algo_node("ANTIKT");
            truthReco->set_input_node("TRUTH");
            
            // Same deal for truth-jet reco: keep JetReco quiet.
            truthReco->Verbosity(0);
            
            se->registerSubsystem(truthReco);
            
            if (vlevel > 0)
                std::cout << "[INFO] (isSim) truth jets: produced node " << truthNode
                << " (R=" << R << ")\n";
        }
    }
    else if (isSim && useDSTTruthJets && vlevel > 0)
    {
        std::cout << "[INFO] (isSim) truth jets: using nodes from DST_JETS (no TruthJetInput reco)"
                  << (auauCandidateSkimOnly ? " [candidate-skim mode]" : "")
                  << "\n";
    }
    
    
    // --------------------------------------------------------------------
    // 5.  Run-information helper (optional but handy)
    // --------------------------------------------------------------------
    if (!isSim)
    {
        auto* trigInfo = new TriggerRunInfoReco();
        trigInfo->Verbosity(vlevel);
        se->registerSubsystem(trigInfo);
    }
    else
    {
        if (vlevel > 0) std::cout << "[isSim] skipping TriggerRunInfoReco" << std::endl;
    }
    
    // Build photon clusters
    class NativeCEMCUESubtractor final : public SubsysReco
    {
    public:
        NativeCEMCUESubtractor(const std::string& name,
                               const std::string& inputNode,
                               const std::string& outputNode)
        : SubsysReco(name)
        , m_inputNode(inputNode)
        , m_outputNode(outputNode)
        {
        }
        
        void setDoSeedExclusion(bool v) { m_doSeedExclusion = v; }
        void setSeedJetNode(const std::string& n) { m_seedJetNode = n; }
        void setSeedMinPt(float v) { m_seedMinPt = v; }
        void setExclusionDR(float v) { m_exclusionDR = v; }
        
        int InitRun(PHCompositeNode* topNode) override
        {
            auto* src = findNode::getClass<TowerInfoContainer>(topNode, m_inputNode);
            if (!src)
            {
                std::cerr << Name() << ": missing input tower node '" << m_inputNode << "'" << std::endl;
                return Fun4AllReturnCodes::ABORTRUN;
            }
            
            PHNodeIterator iter(topNode);
            auto* cemcNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "CEMC"));
            if (!cemcNode)
            {
                std::cerr << Name() << ": missing CEMC composite node" << std::endl;
                return Fun4AllReturnCodes::ABORTRUN;
            }
            
            m_output = findNode::getClass<TowerInfoContainer>(topNode, m_outputNode);
            if (!m_output)
            {
                m_output = dynamic_cast<TowerInfoContainer*>(src->CloneMe());
                if (!m_output)
                {
                    std::cerr << Name() << ": failed to clone input tower container '" << m_inputNode << "'" << std::endl;
                    return Fun4AllReturnCodes::ABORTRUN;
                }
                
                auto* outNode = new PHIODataNode<PHObject>(m_output, m_outputNode, "PHObject");
                cemcNode->addNode(outNode);
            }
            
            return Fun4AllReturnCodes::EVENT_OK;
        }
        
        int process_event(PHCompositeNode* topNode) override
        {
            auto* src = findNode::getClass<TowerInfoContainer>(topNode, m_inputNode);
            auto* dst = findNode::getClass<TowerInfoContainer>(topNode, m_outputNode);
            
            if (!src || !dst)
            {
                std::cerr << Name()
                << ": missing required nodes (src=" << m_inputNode
                << ", dst=" << m_outputNode
                << ")"
                << std::endl;
                return Fun4AllReturnCodes::ABORTEVENT;
            }
            
            ++m_evt;
            
            // --- variantB: build seed-exclusion mask from coresoftware's refined seed jets ---
            std::vector<std::pair<float,float>> seedPositions;  // (eta, phi)
            RawTowerGeomContainer* geomCEMC_excl = nullptr;
            float exclusionDR2 = m_exclusionDR * m_exclusionDR;
            if (m_doSeedExclusion)
            {
                geomCEMC_excl = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
                auto* seeds = findNode::getClass<JetContainer>(topNode, m_seedJetNode);
                if (seeds && geomCEMC_excl)
                {
                    for (auto* jet : *seeds)
                    {
                        if (!jet) continue;
                        if (jet->get_pt() < m_seedMinPt) continue;
                        seedPositions.emplace_back(jet->get_eta(), jet->get_phi());
                    }
                }
                if (Verbosity() > 0)
                {
                    std::cout << "[" << Name() << "] evt=" << m_evt
                    << " seedExclusion: " << seedPositions.size()
                    << " seeds above " << m_seedMinPt << " GeV from " << m_seedJetNode
                    << " (DR=" << m_exclusionDR << ")" << std::endl;
                }
            }
            
            std::vector<double> stripSum;
            std::vector<unsigned int> stripCount;
            std::vector<double> stripMeanCache;
            
            double srcSumEAll = 0.0;
            double srcSumEGood = 0.0;
            unsigned int nSrcFiniteAll = 0;
            unsigned int nSrcFiniteGood = 0;
            unsigned int nSrcPositiveGood = 0;
            
            const unsigned int nchannels = src->size();
            for (unsigned int channel = 0; channel < nchannels; ++channel)
            {
                TowerInfo* srcTower = src->get_tower_at_channel(channel);
                if (!srcTower)
                {
                    continue;
                }
                
                const float energy = srcTower->get_energy();
                if (std::isfinite(energy))
                {
                    srcSumEAll += energy;
                    ++nSrcFiniteAll;
                }
                
                const unsigned int towerkey = src->encode_key(channel);
                const int ieta = src->getTowerEtaBin(towerkey);
                if (ieta < 0)
                {
                    continue;
                }
                
                if (static_cast<std::size_t>(ieta + 1) > stripSum.size())
                {
                    stripSum.resize(static_cast<std::size_t>(ieta + 1), 0.0);
                    stripCount.resize(static_cast<std::size_t>(ieta + 1), 0U);
                }
                
                if (!srcTower->get_isGood())
                {
                    continue;
                }
                
                if (!std::isfinite(energy))
                {
                    continue;
                }
                
                srcSumEGood += energy;
                ++nSrcFiniteGood;
                if (energy > 0.0f) ++nSrcPositiveGood;
                
                // variantB: skip towers within exclusion cone of any seed jet
                if (m_doSeedExclusion && !seedPositions.empty() && geomCEMC_excl)
                {
                    const int iphi = src->getTowerPhiBin(towerkey);
                    const RawTowerDefs::keytype gkey = RawTowerDefs::encode_towerid(
                                                                                    RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);
                    RawTowerGeom* tg = geomCEMC_excl->get_tower_geometry(gkey);
                    if (tg)
                    {
                        const float tEta = tg->get_eta();
                        const float tPhi = tg->get_phi();
                        bool masked = false;
                        for (const auto& sp : seedPositions)
                        {
                            float deta = tEta - sp.first;
                            float dphi = tPhi - sp.second;
                            while (dphi >  M_PI) dphi -= 2.0f * M_PI;
                            while (dphi < -M_PI) dphi += 2.0f * M_PI;
                            if (deta * deta + dphi * dphi < exclusionDR2) { masked = true; break; }
                        }
                        if (masked) continue;
                    }
                }
                
                stripSum.at(static_cast<std::size_t>(ieta)) += energy;
                stripCount.at(static_cast<std::size_t>(ieta)) += 1U;
            }
            
            stripMeanCache.resize(stripSum.size(), 0.0);
            unsigned int nNonEmptyStrips = 0;
            double meanStripMean = 0.0;
            double maxStripMean = -std::numeric_limits<double>::infinity();
            int maxStripEta = -1;
            
            for (std::size_t ieta = 0; ieta < stripSum.size(); ++ieta)
            {
                const unsigned int nstrip = stripCount.at(ieta);
                if (nstrip == 0U) continue;
                
                const double stripMean = stripSum.at(ieta) / static_cast<double>(nstrip);
                stripMeanCache.at(ieta) = stripMean;
                meanStripMean += stripMean;
                ++nNonEmptyStrips;
                
                if (stripMean > maxStripMean)
                {
                    maxStripMean = stripMean;
                    maxStripEta = static_cast<int>(ieta);
                }
            }
            
            if (nNonEmptyStrips > 0)
            {
                meanStripMean /= static_cast<double>(nNonEmptyStrips);
            }
            else
            {
                maxStripMean = 0.0;
            }
            
            double dstSumEAll = 0.0;
            double dstSumEGood = 0.0;
            unsigned int nDstFiniteAll = 0;
            unsigned int nDstFiniteGood = 0;
            unsigned int nDstPositiveGood = 0;
            unsigned int nNegativeGood = 0;
            unsigned int nZeroedGood = 0;
            double mostNegative = 0.0;
            double maxSubtracted = -std::numeric_limits<double>::infinity();
            int maxSubtractedEta = -1;
            float maxSubtractedSrcE = 0.0f;
            float maxSubtractedDstE = 0.0f;
            
            for (unsigned int channel = 0; channel < nchannels; ++channel)
            {
                TowerInfo* srcTower = src->get_tower_at_channel(channel);
                TowerInfo* dstTower = dst->get_tower_at_channel(channel);
                if (!srcTower || !dstTower)
                {
                    continue;
                }
                
                const unsigned int towerkey = src->encode_key(channel);
                const int ieta = src->getTowerEtaBin(towerkey);
                
                float new_energy = 0.0f;
                float stripMean = 0.0f;
                const float src_energy = srcTower->get_energy();
                
                if (ieta >= 0 && static_cast<std::size_t>(ieta) < stripMeanCache.size() && srcTower->get_isGood())
                {
                    stripMean = static_cast<float>(stripMeanCache.at(static_cast<std::size_t>(ieta)));
                    if (std::isfinite(src_energy))
                    {
                        new_energy = src_energy - stripMean;
                    }
                }
                
                dstTower->set_time(srcTower->get_time());
                dstTower->set_energy(new_energy);
                
                if (std::isfinite(new_energy))
                {
                    dstSumEAll += new_energy;
                    ++nDstFiniteAll;
                }
                
                if (srcTower->get_isGood() && std::isfinite(new_energy))
                {
                    dstSumEGood += new_energy;
                    ++nDstFiniteGood;
                    if (new_energy > 0.0f) ++nDstPositiveGood;
                    if (new_energy < 0.0f)
                    {
                        ++nNegativeGood;
                        if (new_energy < mostNegative) mostNegative = new_energy;
                    }
                    if (std::fabs(new_energy) < 1e-6f) ++nZeroedGood;
                }
                
                if (srcTower->get_isGood() && std::isfinite(src_energy))
                {
                    const double subtracted = static_cast<double>(src_energy) - static_cast<double>(new_energy);
                    if (subtracted > maxSubtracted)
                    {
                        maxSubtracted = subtracted;
                        maxSubtractedEta = ieta;
                        maxSubtractedSrcE = src_energy;
                        maxSubtractedDstE = new_energy;
                    }
                }
            }
            
            if (Verbosity() > 0)
            {
                std::cout << "[" << Name() << "] evt=" << m_evt
                << " PHOSUB summary"
                << " | srcNode=" << m_inputNode
                << " dstNode=" << m_outputNode
                << " | nchannels=" << nchannels
                << " | strips(nonEmpty)=" << nNonEmptyStrips
                << " meanStrip=" << std::fixed << std::setprecision(4) << meanStripMean
                << " maxStrip=" << maxStripMean << "@ieta=" << maxStripEta
                << " | srcSumE(all/good)=" << std::setprecision(3) << srcSumEAll << "/" << srcSumEGood
                << " | dstSumE(all/good)=" << dstSumEAll << "/" << dstSumEGood
                << " | goodFinite(src/dst)=" << nSrcFiniteGood << "/" << nDstFiniteGood
                << " | goodE>0(src/dst)=" << nSrcPositiveGood << "/" << nDstPositiveGood
                << " | negGood=" << nNegativeGood
                << " zeroGood=" << nZeroedGood
                << " mostNegative=" << mostNegative
                << " | maxSubtracted=" << maxSubtracted
                << " (src=" << maxSubtractedSrcE
                << " -> dst=" << maxSubtractedDstE
                << ", ieta=" << maxSubtractedEta << ")"
                << std::endl;
            }
            
            if (Verbosity() > 1)
            {
                std::cout << "[" << Name() << "] evt=" << m_evt << " strip means:";
                int printed = 0;
                for (std::size_t ieta = 0; ieta < stripMeanCache.size(); ++ieta)
                {
                    if (stripCount.at(ieta) == 0U) continue;
                    std::cout << " (" << ieta
                    << ": n=" << stripCount.at(ieta)
                    << ", mean=" << std::fixed << std::setprecision(4) << stripMeanCache.at(ieta)
                    << ")";
                    ++printed;
                    if (printed >= 24)
                    {
                        std::cout << " ...";
                        break;
                    }
                }
                std::cout << std::endl;
            }
            
            return Fun4AllReturnCodes::EVENT_OK;
        }
        
    private:
        std::string m_inputNode;
        std::string m_outputNode;
        TowerInfoContainer* m_output{nullptr};
        int m_evt{0};
        bool m_doSeedExclusion{false};
        std::string m_seedJetNode{"AntiKt_TowerInfo_HIRecoSeedsSub_r02"};
        float m_seedMinPt{5.0f};
        float m_exclusionDR{0.4f};
    };
    
    class TowerInfoCanonicalRebaser final : public SubsysReco
    {
    public:
        TowerInfoCanonicalRebaser(const std::string& name,
                                  const std::string& sourceNode,
                                  const std::string& destNode)
        : SubsysReco(name)
        , m_sourceNode(sourceNode)
        , m_destNode(destNode)
        {
        }
        
        int process_event(PHCompositeNode* topNode) override
        {
            auto* src = findNode::getClass<TowerInfoContainer>(topNode, m_sourceNode);
            auto* dst = findNode::getClass<TowerInfoContainer>(topNode, m_destNode);
            
            if (!src || !dst)
            {
                std::cerr << Name()
                << ": missing source/destination tower nodes (src=" << m_sourceNode
                << ", dst=" << m_destNode << ")"
                << std::endl;
                return Fun4AllReturnCodes::ABORTEVENT;
            }
            
            const unsigned int nchannels = std::min(src->size(), dst->size());
            for (unsigned int channel = 0; channel < nchannels; ++channel)
            {
                TowerInfo* srcTower = src->get_tower_at_channel(channel);
                TowerInfo* dstTower = dst->get_tower_at_channel(channel);
                if (!srcTower || !dstTower)
                {
                    continue;
                }
                
                dstTower->set_time(srcTower->get_time());
                dstTower->set_energy(srcTower->get_energy());
            }
            
            return Fun4AllReturnCodes::EVENT_OK;
        }
        
    private:
        std::string m_sourceNode;
        std::string m_destNode;
    };
    
    std::string photonInputClusterNode = "CLUSTERINFO_CEMC";
    bool photonBuilderIsAuAu = false;
    
    if (cfg.clusterUEpipeline == "variantA" && isAuAuLike)
    {
        const std::string nativeCemcNode = towerPrefixPCB + "_CEMC_PHOSUB";
        
        int nativeUEV = 0;
        if (const char* env = std::getenv("RJ_HIUE_VERBOSITY")) nativeUEV = std::atoi(env);
        
        auto* nativeSub = new NativeCEMCUESubtractor("NativeCEMCUESubtractor",
                                                     towerPrefixPCB + "_CEMC",
                                                     nativeCemcNode);
        nativeSub->Verbosity(nativeUEV);
        se->registerSubsystem(nativeSub);
        
        if (nativeUEV > 0)
        {
            auto* auditBefore = new TowerAudit("TowerAudit_PHOSUB_before",
                                               towerPrefixPCB + "_CEMC",
                                               towerPrefixPCB + "_HCALIN_SUB1",
                                               towerPrefixPCB + "_HCALOUT_SUB1",
                                               10);
            auditBefore->Verbosity(nativeUEV);
            se->registerSubsystem(auditBefore);
            
            auto* auditAfter = new TowerAudit("TowerAudit_PHOSUB_after",
                                              nativeCemcNode,
                                              towerPrefixPCB + "_HCALIN_SUB1",
                                              towerPrefixPCB + "_HCALOUT_SUB1",
                                              10);
            auditAfter->Verbosity(nativeUEV);
            se->registerSubsystem(auditAfter);
        }
        
        // Recluster from UE-subtracted PHOSUB towers (ATLAS-like approach:
        // cluster on subtracted input so cluster energies, positions, and
        // tower membership all reflect the subtracted state)
        {
            auto* phoSubClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate_PHOSUB");
            phoSubClusterBuilder->Detector("CEMC");
            phoSubClusterBuilder->set_threshold_energy(0.070);
            const char* _calibroot = std::getenv("CALIBRATIONROOT");
            if (_calibroot && std::string(_calibroot).size())
            {
                std::string emc_prof_phosub = std::string(_calibroot) + "/EmcProfile/CEMCprof_Thresh30MeV.root";
                phoSubClusterBuilder->LoadProfile(emc_prof_phosub);
            }
            phoSubClusterBuilder->set_UseTowerInfo(1);
            phoSubClusterBuilder->set_UseAltZVertex(1);
            phoSubClusterBuilder->setInputTowerNodeName(nativeCemcNode);
            phoSubClusterBuilder->setOutputClusterNodeName("CLUSTERINFO_CEMC_PHOSUB");
            phoSubClusterBuilder->Verbosity(nativeUEV);
            se->registerSubsystem(phoSubClusterBuilder);
        }
        
        photonInputClusterNode = "CLUSTERINFO_CEMC_PHOSUB";
        
        if (vlevel > 0)
        {
            std::cout << "[clusterUEpipeline=variantA] enabled for AuAu-like running"
            << " | nativeCemcNode=" << nativeCemcNode
            << " | photonInputClusterNode=" << photonInputClusterNode
            << " | PhotonClusterBuilder will read PHOSUB clusters + PHOSUB/SUB1 tower nodes"
            << std::endl;
        }
    }
    else if (cfg.clusterUEpipeline == "variantB" && isAuAuLike)
    {
        const std::string nativeCemcNode = towerPrefixPCB + "_CEMC_PHOSUB";
        
        int nativeUEV = 0;
        if (const char* env = std::getenv("RJ_HIUE_VERBOSITY")) nativeUEV = std::atoi(env);
        
        auto* nativeSub = new NativeCEMCUESubtractor("NativeCEMCUESubtractor_B",
                                                     towerPrefixPCB + "_CEMC",
                                                     nativeCemcNode);
        nativeSub->setDoSeedExclusion(true);
        nativeSub->setSeedJetNode("AntiKt_TowerInfo_HIRecoSeedsSub_r02");
        nativeSub->setSeedMinPt(5.0f);
        nativeSub->setExclusionDR(0.4f);
        nativeSub->Verbosity(nativeUEV);
        se->registerSubsystem(nativeSub);
        
        // Recluster from seed-masked UE-subtracted PHOSUB towers
        {
            auto* phoSubClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate_PHOSUB");
            phoSubClusterBuilder->Detector("CEMC");
            phoSubClusterBuilder->set_threshold_energy(0.070);
            const char* _calibroot = std::getenv("CALIBRATIONROOT");
            if (_calibroot && std::string(_calibroot).size())
            {
                std::string emc_prof_phosub = std::string(_calibroot) + "/EmcProfile/CEMCprof_Thresh30MeV.root";
                phoSubClusterBuilder->LoadProfile(emc_prof_phosub);
            }
            phoSubClusterBuilder->set_UseTowerInfo(1);
            phoSubClusterBuilder->set_UseAltZVertex(1);
            phoSubClusterBuilder->setInputTowerNodeName(nativeCemcNode);
            phoSubClusterBuilder->setOutputClusterNodeName("CLUSTERINFO_CEMC_PHOSUB");
            phoSubClusterBuilder->Verbosity(nativeUEV);
            se->registerSubsystem(phoSubClusterBuilder);
        }
        
        photonInputClusterNode = "CLUSTERINFO_CEMC_PHOSUB";
        
        if (vlevel > 0)
        {
            std::cout << "[clusterUEpipeline=variantB] enabled for AuAu-like running"
            << " | nativeCemcNode=" << nativeCemcNode
            << " | seed exclusion: DR=0.4 from HIRecoSeedsSub_r02 (pt>5 GeV)"
            << " | photonInputClusterNode=" << photonInputClusterNode
            << std::endl;
        }
    }
    else if (cfg.clusterUEpipeline == "baseVariant" && isAuAuLike)
    {
        photonBuilderIsAuAu = true;
        
        if (vlevel > 0)
        {
            std::cout << "[clusterUEpipeline=baseVariant] enabled for AuAu-like running"
            << " | PCB m_is_auau=true → isolation uses RETOWER_SUB1/HCAL_SUB1 internally"
            << " | shower shapes use standard CEMC towers"
            << std::endl;
        }
    }
    
    std::string preselectionPhotonNode = "PHOTONCLUSTER_CEMC";
    std::string tightPhotonNode        = "PHOTONCLUSTER_CEMC";
    
    double minPhotonEt = 5.0;
    if (!cfg.jes3_photon_pt_bins.empty())
    {
        auto itMin = std::min_element(cfg.jes3_photon_pt_bins.begin(), cfg.jes3_photon_pt_bins.end());
        if (itMin != cfg.jes3_photon_pt_bins.end() && std::isfinite(*itMin))
        {
            minPhotonEt = *itMin;
        }
    }
    // Replay-foundation capture is deliberately wider than the reader-facing
    // reporting bins.  Keep this as an explicit, fail-closed runtime contract
    // so a 15 GeV JES3/reporting edge cannot silently erase the lower response
    // guard population before RJPhotonCandidateV1 is written.  Direct and
    // writer canary arms receive the same value and therefore remain a valid
    // scientific-neutrality pair.
    if (const char* raw = std::getenv("RJ_REPLAY_PHOTON_CAPTURE_ET_MIN"))
    {
        try
        {
            const double requested = std::stod(detail::trim(std::string(raw)));
            if (!std::isfinite(requested) || requested < 0.0 || requested >= 15.0)
            {
                detail::bail(
                    "RJ_REPLAY_PHOTON_CAPTURE_ET_MIN must be finite, nonnegative, and below the 15 GeV reporting boundary");
            }
            minPhotonEt = requested;
        }
        catch (const std::exception&)
        {
            detail::bail(
                "RJ_REPLAY_PHOTON_CAPTURE_ET_MIN must parse as a finite numeric threshold");
        }
    }
    
    const bool useSamePhotonBDTScores = true;
    const bool usePPG12PPIsoTowerFloor = !isAuAuLike;
    const bool ppg12PhotonYieldPPSim =
        env_truthy_local("RJ_PPG12_PHOTON_YIELD") &&
        isSim && !isAuAuLike;
    // Canonical CEMC shower-shape contract (all systems and sample types):
    // build the 7x7 grid from the complete good-TowerInfo grid.  The pp
    // core/PPG12 contract accepts TowerInfo cells through get_isGood() only;
    // Au+Au-specific status/calibration is already encoded in that flag, so no
    // second analysis-local chi2/CDB mask is applied to the shower grid. pp
    // retains the PPG12 70 MeV cell floor; Au+Au uses a zero configured cell
    // floor. Data and embedded simulation must share both settings. The RawCluster
    // towermap remains a PhotonClusterBuilder diagnostic API only and is never
    // selected implicitly by the unified production path.
    constexpr float kPPG12PPCEMCShapeTowerMinGeV = 0.070f;
    constexpr float kAuAuCEMCShapeTowerMinGeV = 0.0f;
    std::string cemcShowerShapeDiagnosticVariant = "canonical";
    if (const char* raw = std::getenv("RJ_AUAU_SHOWER_SHAPE_DIAGNOSTIC_VARIANT"))
    {
        cemcShowerShapeDiagnosticVariant = detail::trim(std::string(raw));
        std::transform(cemcShowerShapeDiagnosticVariant.begin(),
                       cemcShowerShapeDiagnosticVariant.end(),
                       cemcShowerShapeDiagnosticVariant.begin(),
                       [](unsigned char c){ return std::tolower(c); });
    }
    if (cemcShowerShapeDiagnosticVariant != "canonical" &&
        cemcShowerShapeDiagnosticVariant != "towerinfo70" &&
        cemcShowerShapeDiagnosticVariant != "historical")
    {
        detail::bail(
            "RJ_AUAU_SHOWER_SHAPE_DIAGNOSTIC_VARIANT must be one of "
            "canonical, towerinfo70, or historical; received \"" +
            cemcShowerShapeDiagnosticVariant + "\"");
    }
    if (!isAuAuLike && cemcShowerShapeDiagnosticVariant != "canonical")
    {
        detail::bail(
            "RJ_AUAU_SHOWER_SHAPE_DIAGNOSTIC_VARIANT=" +
            cemcShowerShapeDiagnosticVariant +
            " is an AuAu-only diagnostic and cannot modify the pp contract");
    }

    bool useCoreGoodTowerInfoShapes = true;
    bool useRawClusterTowermapForCEMCShapes = false;
    float resolvedCEMCShapeTowerMinGeV =
        isAuAuLike ? kAuAuCEMCShapeTowerMinGeV : kPPG12PPCEMCShapeTowerMinGeV;
    std::string resolvedCEMCShapeEnergySource = "towerinfo_full_good_grid";
    std::string resolvedCEMCShapeTowerAcceptance = "towerinfo_get_isgood";

    if (isAuAuLike && cemcShowerShapeDiagnosticVariant == "towerinfo70")
    {
        resolvedCEMCShapeTowerMinGeV = kPPG12PPCEMCShapeTowerMinGeV;
    }
    else if (isAuAuLike && cemcShowerShapeDiagnosticVariant == "historical")
    {
        resolvedCEMCShapeTowerMinGeV = kPPG12PPCEMCShapeTowerMinGeV;
        // Reproduce the complete pre-repair AuAu route. Data used the full
        // TowerInfo grid followed by the analysis-local chi2/CDB mask, while
        // embedding read only RawCluster-owned cells. The 70 MeV floor was
        // common to both. This is diagnostic-only and cannot become the
        // default because canonical remains the fail-closed value above.
        useCoreGoodTowerInfoShapes = false;
        resolvedCEMCShapeTowerAcceptance =
            "towerinfo_get_isgood_plus_local_chi2_cdb_mask";
        if (isSimEmbedded)
        {
            useRawClusterTowermapForCEMCShapes = true;
            resolvedCEMCShapeEnergySource = "raw_cluster_towermap_diagnostic";
            resolvedCEMCShapeTowerAcceptance = "raw_cluster_towermap_membership";
        }
    }
    const bool requestedPPG12TruthVertexForReco =
        ppg12PhotonYieldPPSim &&
        (env_truthy_local("RJ_PPG12_PHOTON_YIELD_TRUTH_VERTEX") ||
         env_truthy_local("RJ_PPG12_PHOTON_YIELD_BUILDER_TRUTH_VERTEX") ||
         env_truthy_local("RJ_PPG12_PHOTON_YIELD_RECO_TRUTH_VERTEX"));
    if (requestedPPG12TruthVertexForReco)
    {
        detail::bail(
            "PPG12 truth-for-reconstructed-object vertex mode is forbidden: "
            "reconstructed cluster kinematics/Eiso use reconstructed MBD z; "
            "truth vertices are reserved for SI/DI event weights.");
    }
    const double photonBuilderVzCutCm = std::max(60.0, cfg.vz_cut_cm); // stored reconstruction support; final cuts unchanged
    constexpr float kPPG12PPIsoTowerMin = 0.12f;
    const float photonBuilderIsoTowerMin = usePPG12PPIsoTowerFloor
        ? kPPG12PPIsoTowerMin
        : 0.0f;
    const double recoilJetsIsoTowerMin = isAuAuLike ? 0.0 : cfg.isoTowMin;
    const bool useAuAuTopoClusterIsoCalibration =
        isAuAuLike &&
        env_truthy_local("RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION");
    // The replay writer requires the same p+p topocluster population used by
    // the nominal R=0.4 isolation contract.  Enable it for writer jobs and for
    // their writer-disabled direct controls; otherwise fetchNodes() correctly
    // fails closed on a missing TOPOCLUSTER_ALLCALO node before any candidate,
    // jet, truth, or response rows can be retained.
    const bool useReplayFoundationPPTopoIso =
        !isAuAuLike &&
        (env_truthy_local("RJ_REPLAY_FOUNDATION_V1") ||
         env_truthy_local("RJ_REPLAY_FOUNDATION_CANARY"));
    const bool usePPG12PhotonYieldTopoIso =
        (env_truthy_local("RJ_PPG12_PHOTON_YIELD") &&
         !isAuAuLike &&
         env_bool_local("RJ_PPG12_PHOTON_YIELD_TOPO_ISO", true)) ||
        useAuAuTopoClusterIsoCalibration ||
        useReplayFoundationPPTopoIso;
    const bool ppg12ExcludeCandidateTopo =
        (ppg12PhotonYieldPPSim &&
         env_bool_local("RJ_PPG12_PHOTON_YIELD_EXCLUDE_CANDIDATE_TOPO", false)) ||
        (useAuAuTopoClusterIsoCalibration &&
         env_bool_local("RJ_AUAU_TOPOCLUSTER_EXCLUDE_CANDIDATE", false));

    bool ppg12TopoBuilderRegistered = false;
    auto registerPPG12PhotonYieldTopoBuilder = [&](const char* placement)
    {
        if (ppg12TopoBuilderRegistered)
        {
            return;
        }

        auto* ppg12TopoBuilder = new RawClusterBuilderTopo("RawClusterBuilderTopo_PPG12PhotonYield");
        ppg12TopoBuilder->set_nodename("TOPOCLUSTER_ALLCALO");
        ppg12TopoBuilder->setInputTowerNodePrefix(towerPrefixPCB);
        ppg12TopoBuilder->set_enable_HCal(true);
        ppg12TopoBuilder->set_enable_EMCal(true);
        ppg12TopoBuilder->set_noise(0.0053, 0.0351, 0.0684);
        ppg12TopoBuilder->set_significance(4.0, 2.0, 1.0);
        ppg12TopoBuilder->allow_corner_neighbor(true);
        ppg12TopoBuilder->set_do_split(true);
        ppg12TopoBuilder->set_minE_local_max(1.0, 2.0, 0.5);
        ppg12TopoBuilder->set_R_shower(0.025);
        ppg12TopoBuilder->set_use_only_good_towers(true);
        ppg12TopoBuilder->set_absE(true);
        ppg12TopoBuilder->Verbosity(0);
        se->registerSubsystem(ppg12TopoBuilder);
        ppg12TopoBuilderRegistered = true;

        if (vlevel > 0)
        {
            std::cout << "[PPG12_PHOTON_YIELD_V1] registered pp topo cluster builder"
                      << " node=TOPOCLUSTER_ALLCALO"
                      << " towerPrefix=" << towerPrefixPCB
                      << " ppg12TopoConfig=noise(0.0053,0.0351,0.0684),sig(4,2,1),split,goodTowers,absE"
                      << " placement=" << placement
                      << "\n";
        }
    };

    if (usePPG12PhotonYieldTopoIso)
    {
        registerPPG12PhotonYieldTopoBuilder("before PhotonClusterBuilder");
    }

    auto configurePhotonBuilder =
    [&](PhotonClusterBuilder* builder, const std::string& outNode)
    {
        builder->set_input_cluster_node(photonInputClusterNode);
        builder->set_output_photon_node(outNode);
        builder->set_ET_threshold(static_cast<float>(minPhotonEt));
        builder->set_shower_shape_min_tower_energy(resolvedCEMCShapeTowerMinGeV);
        builder->set_iso_min_tower_energy(photonBuilderIsoTowerMin);
        builder->set_use_ppg12_pp_iso_axis(useSamePhotonBDTScores);
        builder->set_use_ppg12_pp_sim_truth_vertex(false);
        builder->set_use_ppg12_pp_sim_global_mbd_vertex(ppg12PhotonYieldPPSim);
        // CaloAna24/core PhotonClusterBuilder forms the 7x7 shower-shape inputs
        // from the complete good-TowerInfo grid.  Apply that energy-source
        // contract to pp and Au+Au, data and simulation.  In particular, do
        // not let embedded Au+Au silently fall back to the RawCluster towermap:
        // that would give training/response samples different BDT inputs from
        // data reconstructed with the same nominal selection.
        builder->set_use_ppg12_pp_sim_towerinfo_shapes(useCoreGoodTowerInfoShapes);
        builder->set_use_ppg12_topocluster_isolation(usePPG12PhotonYieldTopoIso);
        builder->set_ppg12_topocluster_node("TOPOCLUSTER_ALLCALO");
        builder->set_ppg12_topocluster_iso_radius(0.4f);
        builder->set_ppg12_topocluster_exclude_candidate(ppg12ExcludeCandidateTopo);
        builder->set_skip_ppg12_edge_clusters(useSamePhotonBDTScores);
        builder->set_enable_ss_3x3_moments(isAuAuLike);
        builder->set_use_raw_cluster_towermap_for_cemc_shapes(
            useRawClusterTowermapForCEMCShapes);
        
        builder->set_use_vz_cut(cfg.use_vz_cut);
        builder->set_vz_cut_cm(photonBuilderVzCutCm);
        
        builder->set_is_auau(photonBuilderIsAuAu);
        if ((cfg.clusterUEpipeline == "variantA" || cfg.clusterUEpipeline == "variantB") && isAuAuLike)
        {
            builder->set_emc_tower_node(towerPrefixPCB + "_CEMC_PHOSUB");
            builder->set_ihcal_tower_node(towerPrefixPCB + "_HCALIN_SUB1");
            builder->set_ohcal_tower_node(towerPrefixPCB + "_HCALOUT_SUB1");
        }
        else if (cfg.clusterUEpipeline == "baseVariant" && isAuAuLike)
        {
            builder->set_emc_tower_node(towerPrefixPCB + "_CEMC");
            builder->set_ihcal_tower_node(towerPrefixPCB + "_HCALIN");
            builder->set_ohcal_tower_node(towerPrefixPCB + "_HCALOUT");
            builder->set_tower_node_prefix(towerPrefixPCB);
        }
        else if (isAuAuLike)
        {
            // noSub: standard towers, no AuAu iso path
            builder->set_emc_tower_node(towerPrefixPCB + "_CEMC");
            builder->set_ihcal_tower_node(towerPrefixPCB + "_HCALIN");
            builder->set_ohcal_tower_node(towerPrefixPCB + "_HCALOUT");
        }
        
        builder->Verbosity(vlevel);
    };
    
    auto* photonBuilder = new PhotonClusterBuilder("PhotonClusterBuilder");
    configurePhotonBuilder(photonBuilder, "PHOTONCLUSTER_CEMC");
    if (!useSamePhotonBDTScores)
    {
        se->registerSubsystem(photonBuilder);
    }
    
    if (attachPPNPBScore)
    {
        if (cfg.npb_model_file.empty())
        {
            detail::bail("NPB preselection variants require npb_model_file in analysis_config.yaml");
        }
        if (cfg.npb_features.empty())
        {
            detail::bail("NPB preselection variants require npb_features in analysis_config.yaml");
        }

        if (useSamePhotonBDTScores)
        {
            float ppNPBScoreMinEt = 6.0f;
            if (const char* env = std::getenv("RJ_PP_NPB_SCORE_MIN_ET"))
            {
                char* end = nullptr;
                const float parsed = std::strtof(env, &end);
                if (end != env && std::isfinite(parsed)) ppNPBScoreMinEt = parsed;
            }
            float ppNPBScoreMaxEt = 40.0f;
            if (const char* env = std::getenv("RJ_PP_NPB_SCORE_MAX_ET"))
            {
                char* end = nullptr;
                const float parsed = std::strtof(env, &end);
                if (end != env && (std::isfinite(parsed) || std::isinf(parsed))) ppNPBScoreMaxEt = parsed;
            }
            preselectionPhotonNode = "PHOTONCLUSTER_CEMC";
            photonBuilder->add_named_bdt_score("npb_score",
                                               cfg.npb_model_file,
                                               cfg.npb_features,
                                               ppNPBScoreMinEt,
                                               ppNPBScoreMaxEt,
                                               0.7f);
            if (ppPhotonIDTrainingWantsNPBAudit && !fanoutUsesNPB)
            {
                std::cout << "[PPPhotonIDTrainingTree] attaching pp NPB audit score on PHOTONCLUSTER_CEMC"
                          << " while keeping preselection=" << cfg.preselection
                          << " (score stored for downstream nbkg_cut validation, not used as the preselection gate)\n";
            }
        }
        else
        {
            preselectionPhotonNode = "PHOTONCLUSTER_CEMC_NPB";

            auto* photonBuilderNPB = new PhotonClusterBuilder("PhotonClusterBuilder_NPB");
            configurePhotonBuilder(photonBuilderNPB, preselectionPhotonNode);
            photonBuilderNPB->set_do_bdt(true);
            photonBuilderNPB->set_bdt_model_file(cfg.npb_model_file);
            photonBuilderNPB->set_bdt_feature_list(cfg.npb_features);
            se->registerSubsystem(photonBuilderNPB);
        }
    }

    if (fanoutUsesAuAuNPB)
    {
        if (cfg.auau_npb_model_file.empty())
        {
            detail::bail("preselection=auauOnlyNPB requires auau_npb_model_file in analysis_config.yaml");
        }
        if (cfg.auau_npb_features.empty())
        {
            detail::bail("preselection=auauOnlyNPB requires auau_npb_features in analysis_config.yaml");
        }

        preselectionPhotonNode = "PHOTONCLUSTER_CEMC";
        photonBuilder->add_named_bdt_score("auau_npb_score",
                                           cfg.auau_npb_model_file,
                                           cfg.auau_npb_features,
                                           5.0f,
                                           80.0f,
                                           static_cast<float>(cfg.photon_eta_abs_max));
    }
    
    if (fanoutUsesNewPPG12Tight)
    {
        if (cfg.tight_bdt_model_file.empty())
        {
            detail::bail("tight=newPPG12 requires tight_bdt_model_file in analysis_config.yaml");
        }
        if (cfg.tight_bdt_features.empty())
        {
            detail::bail("tight=newPPG12 requires tight_bdt_features in analysis_config.yaml");
        }

        if (useSamePhotonBDTScores)
        {
            tightPhotonNode = "PHOTONCLUSTER_CEMC";
            photonBuilder->add_named_bdt_score("tight_bdt_score",
                                               cfg.tight_bdt_model_file,
                                               cfg.tight_bdt_features,
                                               7.0f);
            // Deployed PPG12 keeps both score branches in the SlimTree, then
            // uses the once-smeared reconstructed ET only to choose the branch:
            // base_v3E for 8 <= ET < 35 GeV and base_E otherwise.  Evaluate
            // both models here from the original reconstructed feature row;
            // RecoilJets performs the route after materializing that ET.
            if (ppg12PhotonYieldPPSim)
            {
                std::string baseEModelFile = cfg.ppg12_base_e_model_file;
                if (baseEModelFile.empty())
                {
                    baseEModelFile = cfg.tight_bdt_model_file;
                    const std::string baseV3EToken = "base_v3E";
                    const std::size_t modelTokenPos = baseEModelFile.find(baseV3EToken);
                    if (modelTokenPos == std::string::npos)
                    {
                        detail::bail(
                            "PPG12 pp-SIM photon-yield mode requires either "
                            "ppg12_base_e_model_file or a base_v3E "
                            "tight_bdt_model_file whose legacy sibling path "
                            "can be resolved; received " + baseEModelFile);
                    }
                    baseEModelFile.replace(
                        modelTokenPos, baseV3EToken.size(), "base_E");
                }
                const std::vector<std::string> baseEFeatures = {
                    "cluster_Et", "vertexz", "cluster_Eta", "e11_over_e33",
                    "cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4"};
                photonBuilder->add_named_bdt_score("tight_bdt_score_base_e",
                                                   baseEModelFile,
                                                   baseEFeatures);
            }

            // The replay foundation keeps the campaign classifier and the
            // historical PPG12 classifier as distinct model evaluations on
            // the identical loose candidate.  The reference model is never
            // allowed to drive the nominal tight selection.
            if (!isAuAuLike)
            {
                const char* referenceModel = std::getenv("RJ_REPLAY_REFERENCE_MODEL_FILE");
                if (referenceModel && *referenceModel)
                {
                    photonBuilder->add_named_bdt_score("ppg12_reference_bdt_score",
                                                       referenceModel,
                                                       cfg.tight_bdt_features,
                                                       7.0f);
                }
            }
        }
        else
        {
            tightPhotonNode = "PHOTONCLUSTER_CEMC_TIGHTBDT";

            auto* photonBuilderTightBDT = new PhotonClusterBuilder("PhotonClusterBuilder_TightBDT");
            configurePhotonBuilder(photonBuilderTightBDT, tightPhotonNode);
            photonBuilderTightBDT->set_do_bdt(true);
            photonBuilderTightBDT->set_bdt_model_file(cfg.tight_bdt_model_file);
            photonBuilderTightBDT->set_bdt_feature_list(cfg.tight_bdt_features);
            se->registerSubsystem(photonBuilderTightBDT);
        }
    }

    auto featureListOrFallback = [](const std::vector<std::string>& primary,
                                    const std::vector<std::string>& fallback) -> std::vector<std::string>
    {
        return primary.empty() ? fallback : primary;
    };
    auto appendCentralityFeature = [](std::vector<std::string> features) -> std::vector<std::string>
    {
        const auto has = std::find(features.begin(), features.end(), "centrality") != features.end() ||
                         std::find(features.begin(), features.end(), "cent") != features.end();
        if (!has) features.push_back("centrality");
        return features;
    };
    // Fail loud when a tight mode's OWN feature key is missing. The exported BDT is
    // positional (feature names stripped to f0..fN at export), so silently falling
    // back to a DIFFERENT feature set corrupts every score with no error. Modes
    // whose fallback would resolve to a different shower-shape set (3x3 / base3x3)
    // must require their own key instead of guessing.
    auto requireFeatures = [](const std::vector<std::string>& primary,
                              const char* keyName,
                              const std::string& mode) -> const std::vector<std::string>&
    {
        if (primary.empty())
        {
            std::ostringstream oss;
            oss << "tightMode '" << mode << "' requires config key '" << keyName
                << "' but it is missing/empty. Refusing to silently fall back to a different "
                   "feature set (the model is positional; a wrong feature list corrupts every "
                   "score). Define '" << keyName << "' in the active analysis config.";
            throw std::runtime_error(oss.str());
        }
        return primary;
    };
    auto validateModelCount = [](const std::string& label, std::size_t got, std::size_t expected)
    {
        if (got == expected) return;
        std::ostringstream msg;
        msg << label << " requires " << expected << " model file(s); got " << got;
        detail::bail(msg.str());
    };
    auto fmtInt3 = [](double x) -> std::string
    {
        std::ostringstream os;
        os << std::setw(3) << std::setfill('0') << static_cast<int>(std::llround(x));
        return os.str();
    };
    auto ptTag = [&](double lo, double hi) -> std::string
    {
        return std::string("pt_") + fmtInt3(lo) + "_" + fmtInt3(hi);
    };
    auto centTag = [&](int lo, int hi) -> std::string
    {
        return std::string("cent_") + fmtInt3(lo) + "_" + fmtInt3(hi);
    };
    auto expandedModelPath = [&](const std::string& modelId) -> std::string
    {
        if (cfg.auau_tight_bdt_expanded_model_dir.empty()) return std::string{};
        return cfg.auau_tight_bdt_expanded_model_dir + "/auau_tight_bdt_" + modelId + "_tmva.root";
    };
    auto mlpModelPath = [&](const std::string& filename) -> std::string
    {
        if (cfg.auau_tight_mlp_model_dir.empty()) return std::string{};
        return cfg.auau_tight_mlp_model_dir + "/" + filename;
    };
    auto expandedCentModels = [&](const std::string& product,
                                  const std::string& suffix,
                                  const std::vector<int>& edges) -> std::vector<std::string>
    {
        std::vector<std::string> files;
        if (cfg.auau_tight_bdt_expanded_model_dir.empty() || edges.size() < 2) return files;
        for (std::size_t i = 0; i + 1 < edges.size(); ++i)
        {
            files.push_back(expandedModelPath(product + "_" + suffix + "_" + centTag(edges[i], edges[i + 1])));
        }
        return files;
    };
    auto expandedPtModels = [&](const std::string& product,
                                const std::vector<double>& ptEdges) -> std::vector<std::string>
    {
        std::vector<std::string> files;
        if (cfg.auau_tight_bdt_expanded_model_dir.empty() || ptEdges.size() < 2) return files;
        for (std::size_t i = 0; i + 1 < ptEdges.size(); ++i)
        {
            files.push_back(expandedModelPath(product + "_" + ptTag(ptEdges[i], ptEdges[i + 1])));
        }
        return files;
    };
    auto expandedPtCentModels = [&](const std::string& product,
                                    const std::vector<double>& ptEdges,
                                    const std::vector<int>& centEdges) -> std::vector<std::string>
    {
        std::vector<std::string> files;
        if (cfg.auau_tight_bdt_expanded_model_dir.empty() || ptEdges.size() < 2 || centEdges.size() < 2) return files;
        for (std::size_t ip = 0; ip + 1 < ptEdges.size(); ++ip)
        {
            for (std::size_t ic = 0; ic + 1 < centEdges.size(); ++ic)
            {
                files.push_back(expandedModelPath(product + "_" + ptTag(ptEdges[ip], ptEdges[ip + 1]) + "_" + centTag(centEdges[ic], centEdges[ic + 1])));
            }
        }
        return files;
    };

    struct AuAuBDTRuntimeConfig
    {
        std::string modelFile;
        std::vector<std::string> features;
        std::vector<int> centEdges;
        std::vector<std::string> centModelFiles;
        std::vector<double> ptEdges;
        std::vector<std::string> ptModelFiles;
        std::vector<std::string> ptCentModelFiles;
        std::string ptFallbackModelFile;
        std::vector<std::string> ptFallbackCentModelFiles;
        double ptFallbackMin = 35.0;
        double ptFallbackMax = 40.0;
        double applyPtMin = std::numeric_limits<double>::quiet_NaN();
        double applyPtMax = std::numeric_limits<double>::quiet_NaN();
        std::vector<std::string> workingPointEntries;
        bool usesBuilderScore = false;
        bool active = false;
    };

    auto resolveAuAuBDTRuntimeConfig = [&](const std::string& tightMode) -> AuAuBDTRuntimeConfig
    {
        AuAuBDTRuntimeConfig out;
        if (!yamlcfg::IsAuAuTightBDTMode(tightMode)) return out;
        out.active = true;
        out.features = cfg.auau_tight_bdt_features;
        out.applyPtMin = cfg.auau_tight_bdt_apply_pt_min;
        out.applyPtMax = cfg.auau_tight_bdt_apply_pt_max;
        out.workingPointEntries = cfg.auau_tight_bdt_working_point_entries;

        if (tightMode == "auauEmbeddedBDT")
        {
            out.modelFile = cfg.auau_tight_bdt_model_file;
            out.usesBuilderScore = true;
        }
        else if (tightMode == "centINDcontrol")
        {
            out.modelFile = cfg.auau_tight_bdt_centINDcontrol_model_file.empty()
                ? cfg.auau_tight_bdt_model_file
                : cfg.auau_tight_bdt_centINDcontrol_model_file;
            out.features = featureListOrFallback(cfg.auau_tight_bdt_centINDcontrol_features, cfg.auau_tight_bdt_features);
        }
        else if (tightMode == "centAsFeat")
        {
            out.modelFile = cfg.auau_tight_bdt_centAsFeat_model_file.empty()
                ? cfg.auau_tight_bdt_model_file
                : cfg.auau_tight_bdt_centAsFeat_model_file;
            out.features = appendCentralityFeature(featureListOrFallback(cfg.auau_tight_bdt_centAsFeat_features, cfg.auau_tight_bdt_features));
        }
        else if (tightMode == "centDepBDTs")
        {
            out.features = featureListOrFallback(cfg.auau_tight_bdt_centDep_features, cfg.auau_tight_bdt_features);
            out.centEdges = cfg.centrality_edges;
            validateModelCount("tight=centDepBDTs", cfg.auau_tight_bdt_centDep_model_files.size(),
                               out.centEdges.size() >= 2 ? out.centEdges.size() - 1 : 0);
            out.centModelFiles = cfg.auau_tight_bdt_centDep_model_files;
        }
        else if (tightMode == "auauNoCentBDT")
        {
            out.modelFile = cfg.auau_tight_bdt_noCent_model_file.empty()
                ? expandedModelPath("centINDcontrol_pt5to40")
                : cfg.auau_tight_bdt_noCent_model_file;
            out.features = featureListOrFallback(cfg.auau_tight_bdt_centINDcontrol_features, cfg.auau_tight_bdt_features);
        }
        else if (tightMode == "auauCentInputBDT")
        {
            out.modelFile = cfg.auau_tight_bdt_centInput_model_file.empty()
                ? expandedModelPath("centAsFeat_pt5to40")
                : cfg.auau_tight_bdt_centInput_model_file;
            out.features = appendCentralityFeature(featureListOrFallback(cfg.auau_tight_bdt_centAsFeat_features, cfg.auau_tight_bdt_features));
        }
        else if (tightMode == "auauCentInput3x3BDT")
        {
            out.modelFile = cfg.auau_tight_bdt_centInput3x3_model_file.empty()
                ? expandedModelPath("centAsFeat3x3_pt5to40")
                : cfg.auau_tight_bdt_centInput3x3_model_file;
            out.features = appendCentralityFeature(requireFeatures(cfg.auau_tight_bdt_centAsFeat3x3_features, "auau_tight_bdt_centAsFeat3x3_features", tightMode));
        }
        else if (tightMode == "auauCentInputBase3x3BDT")
        {
            out.modelFile = cfg.auau_tight_bdt_centInputBase3x3_model_file.empty()
                ? expandedModelPath("centAsFeatBase3x3_pt15to30")
                : cfg.auau_tight_bdt_centInputBase3x3_model_file;
            out.features = appendCentralityFeature(requireFeatures(cfg.auau_tight_bdt_centAsFeatBase3x3_features, "auau_tight_bdt_centAsFeatBase3x3_features", tightMode));
        }
        else if (tightMode == "auauCentInputMinOptBDT")
        {
            out.modelFile = cfg.auau_tight_bdt_centInputMinOpt_model_file.empty()
                ? expandedModelPath("centAsFeatMinOpt_pt5to40")
                : cfg.auau_tight_bdt_centInputMinOpt_model_file;
            out.features = appendCentralityFeature(featureListOrFallback(cfg.auau_tight_bdt_centAsFeat_features, cfg.auau_tight_bdt_features));
        }
        else if (tightMode == "auauCent3BDT" || tightMode == "auauCent7BDT")
        {
            out.features = featureListOrFallback(cfg.auau_tight_bdt_centDep_features, cfg.auau_tight_bdt_features);
            out.centEdges = (tightMode == "auauCent3BDT") ? cfg.auau_tight_bdt_cent3_edges : cfg.auau_tight_bdt_cent7_edges;
            out.centModelFiles = (tightMode == "auauCent3BDT") ? cfg.auau_tight_bdt_cent3_model_files : cfg.auau_tight_bdt_cent7_model_files;
            if (out.centModelFiles.empty())
            {
                out.centModelFiles = expandedCentModels(tightMode == "auauCent3BDT" ? "centDepBDTs" : "centDepFineBDTs",
                                                        "pt5to40",
                                                        out.centEdges);
            }
            validateModelCount("tight=" + tightMode, out.centModelFiles.size(),
                               out.centEdges.size() >= 2 ? out.centEdges.size() - 1 : 0);
        }
        else if (tightMode == "auauPtBinCentInputBDT")
        {
            out.features = appendCentralityFeature(featureListOrFallback(cfg.auau_tight_bdt_centAsFeat_features, cfg.auau_tight_bdt_features));
            out.ptEdges = cfg.auau_tight_bdt_pt_bin_edges;
            out.ptModelFiles = cfg.auau_tight_bdt_ptBinCentInput_model_files;
            if (out.ptModelFiles.empty())
            {
                out.ptModelFiles = expandedPtModels("ptBinCentAsFeat", out.ptEdges);
            }
            out.ptFallbackModelFile = cfg.auau_tight_bdt_ptBinCentInput_fallback_model_file.empty()
                ? expandedModelPath("centAsFeat_pt5to40")
                : cfg.auau_tight_bdt_ptBinCentInput_fallback_model_file;
            out.ptFallbackMin = cfg.auau_tight_bdt_pt_fallback_min;
            out.ptFallbackMax = cfg.auau_tight_bdt_pt_fallback_max;
            validateModelCount("tight=auauPtBinCentInputBDT", out.ptModelFiles.size(),
                               out.ptEdges.size() >= 2 ? out.ptEdges.size() - 1 : 0);
        }
        else if (tightMode == "auauPtCent3BDT" || tightMode == "auauPtCent7BDT")
        {
            out.features = featureListOrFallback(cfg.auau_tight_bdt_centDep_features, cfg.auau_tight_bdt_features);
            out.ptEdges = cfg.auau_tight_bdt_pt_bin_edges;
            out.centEdges = (tightMode == "auauPtCent3BDT") ? cfg.auau_tight_bdt_cent3_edges : cfg.auau_tight_bdt_cent7_edges;
            out.ptCentModelFiles = (tightMode == "auauPtCent3BDT") ? cfg.auau_tight_bdt_ptCent3_model_files : cfg.auau_tight_bdt_ptCent7_model_files;
            out.ptFallbackCentModelFiles = (tightMode == "auauPtCent3BDT") ? cfg.auau_tight_bdt_ptCent3_fallback_model_files : cfg.auau_tight_bdt_ptCent7_fallback_model_files;
            if (out.ptCentModelFiles.empty())
            {
                out.ptCentModelFiles = expandedPtCentModels(tightMode == "auauPtCent3BDT" ? "ptCentDep3" : "ptCentDepFine",
                                                            out.ptEdges,
                                                            out.centEdges);
            }
            if (out.ptFallbackCentModelFiles.empty())
            {
                out.ptFallbackCentModelFiles = expandedCentModels(tightMode == "auauPtCent3BDT" ? "centDepBDTs" : "centDepFineBDTs",
                                                                  "pt5to40",
                                                                  out.centEdges);
            }
            out.ptFallbackMin = cfg.auau_tight_bdt_pt_fallback_min;
            out.ptFallbackMax = cfg.auau_tight_bdt_pt_fallback_max;
            const std::size_t nPt = out.ptEdges.size() >= 2 ? out.ptEdges.size() - 1 : 0;
            const std::size_t nCent = out.centEdges.size() >= 2 ? out.centEdges.size() - 1 : 0;
            validateModelCount("tight=" + tightMode, out.ptCentModelFiles.size(), nPt * nCent);
            validateModelCount("tight=" + tightMode + " fallback", out.ptFallbackCentModelFiles.size(), nCent);
        }
        else if (tightMode == "auauEtFineCentInputBDT")
        {
            out.features = appendCentralityFeature(requireFeatures(cfg.auau_tight_bdt_centAsFeatBase3x3_features, "auau_tight_bdt_centAsFeatBase3x3_features", tightMode));
            out.ptEdges = cfg.auau_tight_bdt_etfine_pt_bin_edges.empty() ? cfg.auau_tight_bdt_pt_bin_edges : cfg.auau_tight_bdt_etfine_pt_bin_edges;
            out.ptModelFiles = cfg.auau_tight_bdt_etFineCentInput_model_files;
            if (out.ptModelFiles.empty())
            {
                out.ptModelFiles = expandedPtModels("ptFine_centInput", out.ptEdges);
            }
            validateModelCount("tight=auauEtFineCentInputBDT", out.ptModelFiles.size(),
                               out.ptEdges.size() >= 2 ? out.ptEdges.size() - 1 : 0);
        }
        else if (tightMode == "auauEtFineCent3BDT" || tightMode == "auauEtFineCent7BDT")
        {
            out.features = featureListOrFallback(cfg.auau_tight_bdt_centDep_features, cfg.auau_tight_bdt_features);
            out.ptEdges = cfg.auau_tight_bdt_etfine_pt_bin_edges.empty() ? cfg.auau_tight_bdt_pt_bin_edges : cfg.auau_tight_bdt_etfine_pt_bin_edges;
            out.centEdges = (tightMode == "auauEtFineCent3BDT") ? cfg.auau_tight_bdt_cent3_edges : cfg.auau_tight_bdt_cent7_edges;
            out.ptCentModelFiles = (tightMode == "auauEtFineCent3BDT") ? cfg.auau_tight_bdt_etFineCent3_model_files : cfg.auau_tight_bdt_etFineCent7_model_files;
            if (out.ptCentModelFiles.empty())
            {
                const std::string defaultProduct = tightMode == "auauEtFineCent3BDT" ? "ptFine_cent3" : "ptFine_cent7";
                const std::string configuredProduct = tightMode == "auauEtFineCent3BDT"
                    ? cfg.auau_tight_bdt_etFineCent3_product
                    : cfg.auau_tight_bdt_etFineCent7_product;
                const std::string product = configuredProduct.empty() ? defaultProduct : configuredProduct;
                out.ptCentModelFiles = expandedPtCentModels(product,
                                                            out.ptEdges,
                                                            out.centEdges);
            }
            const std::size_t nPt = out.ptEdges.size() >= 2 ? out.ptEdges.size() - 1 : 0;
            const std::size_t nCent = out.centEdges.size() >= 2 ? out.centEdges.size() - 1 : 0;
            validateModelCount("tight=" + tightMode, out.ptCentModelFiles.size(), nPt * nCent);
        }

        if (out.features.empty())
        {
            detail::bail("AuAu tight BDT mode " + tightMode + " requires a non-empty feature list in analysis_config.yaml");
        }
        if (out.modelFile.empty() && out.centModelFiles.empty() && out.ptModelFiles.empty() && out.ptCentModelFiles.empty())
        {
            detail::bail("AuAu tight BDT mode " + tightMode + " requires model file(s) in analysis_config.yaml");
        }
        if ((tightMode == "auauPtBinCentInputBDT") && out.ptFallbackModelFile.empty())
        {
            detail::bail("tight=auauPtBinCentInputBDT requires auau_tight_bdt_ptBinCentInput_fallback_model_file");
        }
        return out;
    };

    struct AuAuMLPRuntimeConfig
    {
        std::string modelFile;
        double minIntercept = 0.80;
        double minSlope = 0.0;
        double maxScore = 1.0;
        double nonTightMinIntercept = 0.20;
        double nonTightMinSlope = 0.0;
        double nonTightMaxIntercept = 0.80;
        double nonTightMaxSlope = 0.0;
        double applyPtMin = std::numeric_limits<double>::quiet_NaN();
        double applyPtMax = std::numeric_limits<double>::quiet_NaN();
        std::vector<std::string> workingPointEntries;
        bool active = false;
    };

    auto resolveAuAuMLPRuntimeConfig = [&](const std::string& tightMode) -> AuAuMLPRuntimeConfig
    {
        AuAuMLPRuntimeConfig out;
        if (!yamlcfg::IsAuAuTightMLPMode(tightMode)) return out;
        out.active = true;
        out.minIntercept = cfg.auau_tight_mlp_min_intercept;
        out.minSlope = cfg.auau_tight_mlp_min_slope;
        out.maxScore = cfg.auau_tight_mlp_max;
        out.nonTightMinIntercept = cfg.auau_nontight_mlp_min_intercept;
        out.nonTightMinSlope = cfg.auau_nontight_mlp_min_slope;
        out.nonTightMaxIntercept = cfg.auau_nontight_mlp_max_intercept;
        out.nonTightMaxSlope = cfg.auau_nontight_mlp_max_slope;
        out.applyPtMin = cfg.auau_tight_mlp_apply_pt_min;
        out.applyPtMax = cfg.auau_tight_mlp_apply_pt_max;
        out.workingPointEntries = cfg.auau_tight_mlp_working_point_entries;
        if (tightMode == "auauCentInputMLP")
        {
            out.modelFile = cfg.auau_tight_mlp_centInput_model_file.empty()
                ? mlpModelPath("auau_tight_mlp_centInput_pt1535.json")
                : cfg.auau_tight_mlp_centInput_model_file;
        }
        else if (tightMode == "auauNoCentBase3x3MLP")
        {
            out.modelFile = cfg.auau_tight_mlp_noCentBase3x3_model_file.empty()
                ? mlpModelPath("auau_tight_mlp_noCentBase3x3_pt1535.json")
                : cfg.auau_tight_mlp_noCentBase3x3_model_file;
        }
        else if (tightMode == "auauCentInputBase3x3MLP")
        {
            out.modelFile = cfg.auau_tight_mlp_centInputBase3x3_model_file.empty()
                ? mlpModelPath("auau_tight_mlp_centInputBase3x3_pt1535.json")
                : cfg.auau_tight_mlp_centInputBase3x3_model_file;
        }
        else if (tightMode == "auauHighPtDistilledKitchenMLP")
        {
            out.modelFile = cfg.auau_tight_mlp_model_file.empty()
                ? mlpModelPath("auau_tight_mlp_highPtDistilledKitchen_v2.json")
                : cfg.auau_tight_mlp_model_file;
        }
        if (out.modelFile.empty() && !cfg.auau_tight_mlp_model_file.empty())
        {
            out.modelFile = cfg.auau_tight_mlp_model_file;
        }
        if (out.modelFile.empty())
        {
            detail::bail("AuAu tight MLP mode " + tightMode + " requires auau_tight_mlp_model_dir or an explicit model file");
        }
        return out;
    };

    struct AuAuBDTMLPStackRuntimeConfig
    {
        std::string modelFile;
        double minIntercept = 0.80;
        double minSlope = 0.0;
        double maxScore = 1.0;
        double nonTightMinIntercept = 0.20;
        double nonTightMinSlope = 0.0;
        double nonTightMaxIntercept = 0.80;
        double nonTightMaxSlope = 0.0;
        double applyPtMin = std::numeric_limits<double>::quiet_NaN();
        double applyPtMax = std::numeric_limits<double>::quiet_NaN();
        std::vector<std::string> workingPointEntries;
        std::string bdtMode;
        std::string mlpMode;
        bool active = false;
    };

    auto resolveAuAuBDTMLPStackRuntimeConfig = [&](const std::string& tightMode) -> AuAuBDTMLPStackRuntimeConfig
    {
        AuAuBDTMLPStackRuntimeConfig out;
        if (!yamlcfg::IsAuAuTightBDTMLPStackMode(tightMode)) return out;
        out.active = true;
        out.modelFile = cfg.auau_tight_bdt_mlp_stack_model_file;
        out.minIntercept = cfg.auau_tight_bdt_mlp_stack_min_intercept;
        out.minSlope = cfg.auau_tight_bdt_mlp_stack_min_slope;
        out.maxScore = cfg.auau_tight_bdt_mlp_stack_max;
        out.nonTightMinIntercept = cfg.auau_nontight_bdt_mlp_stack_min_intercept;
        out.nonTightMinSlope = cfg.auau_nontight_bdt_mlp_stack_min_slope;
        out.nonTightMaxIntercept = cfg.auau_nontight_bdt_mlp_stack_max_intercept;
        out.nonTightMaxSlope = cfg.auau_nontight_bdt_mlp_stack_max_slope;
        out.applyPtMin = cfg.auau_tight_bdt_mlp_stack_apply_pt_min;
        out.applyPtMax = cfg.auau_tight_bdt_mlp_stack_apply_pt_max;
        out.workingPointEntries = cfg.auau_tight_bdt_mlp_stack_working_point_entries;
        out.bdtMode = cfg.auau_tight_bdt_mlp_stack_bdt_mode.empty() ? "auauEtFineCent7BDT" : cfg.auau_tight_bdt_mlp_stack_bdt_mode;
        out.mlpMode = cfg.auau_tight_bdt_mlp_stack_mlp_mode.empty() ? "auauCentInputBase3x3MLP" : cfg.auau_tight_bdt_mlp_stack_mlp_mode;
        if (out.modelFile.empty())
        {
            detail::bail("AuAu BDT+MLP stack mode requires auau_tight_bdt_mlp_stack_model_file");
        }
        if (!yamlcfg::IsAuAuTightBDTMode(out.bdtMode))
        {
            detail::bail("AuAu BDT+MLP stack bdt_mode must be an AuAu tight BDT mode");
        }
        if (!yamlcfg::IsAuAuTightMLPMode(out.mlpMode))
        {
            detail::bail("AuAu BDT+MLP stack mlp_mode must be an AuAu tight MLP mode");
        }
        return out;
    };

    struct AuAuLogRegRuntimeConfig
    {
        std::string modelFile;
        double minIntercept = 0.80;
        double minSlope = 0.0;
        double maxScore = 1.0;
        double nonTightMinIntercept = 0.20;
        double nonTightMinSlope = 0.0;
        double nonTightMaxIntercept = 0.80;
        double nonTightMaxSlope = 0.0;
        double applyPtMin = std::numeric_limits<double>::quiet_NaN();
        double applyPtMax = std::numeric_limits<double>::quiet_NaN();
        std::vector<std::string> workingPointEntries;
        bool active = false;
    };

    auto resolveAuAuLogRegRuntimeConfig = [&](const std::string& tightMode) -> AuAuLogRegRuntimeConfig
    {
        AuAuLogRegRuntimeConfig out;
        if (!yamlcfg::IsAuAuTightLogRegMode(tightMode)) return out;
        out.active = true;
        out.modelFile = cfg.auau_tight_logreg_model_file;
        out.minIntercept = cfg.auau_tight_logreg_min_intercept;
        out.minSlope = cfg.auau_tight_logreg_min_slope;
        out.maxScore = cfg.auau_tight_logreg_max;
        out.nonTightMinIntercept = cfg.auau_nontight_logreg_min_intercept;
        out.nonTightMinSlope = cfg.auau_nontight_logreg_min_slope;
        out.nonTightMaxIntercept = cfg.auau_nontight_logreg_max_intercept;
        out.nonTightMaxSlope = cfg.auau_nontight_logreg_max_slope;
        out.applyPtMin = cfg.auau_tight_logreg_apply_pt_min;
        out.applyPtMax = cfg.auau_tight_logreg_apply_pt_max;
        out.workingPointEntries = cfg.auau_tight_logreg_working_point_entries;
        if (out.modelFile.empty())
        {
            detail::bail("AuAu tight logistic-regression mode requires auau_tight_logreg_model_file");
        }
        return out;
    };

    std::vector<std::string> auauCentDepScoreNames;
    const AuAuBDTMLPStackRuntimeConfig leaderAuAuStackRuntime = resolveAuAuBDTMLPStackRuntimeConfig(cfg.tight);
    const AuAuLogRegRuntimeConfig leaderAuAuLogRegRuntime = resolveAuAuLogRegRuntimeConfig(cfg.tight);
    const AuAuBDTRuntimeConfig leaderAuAuRuntime = resolveAuAuBDTRuntimeConfig(leaderAuAuStackRuntime.active ? leaderAuAuStackRuntime.bdtMode : cfg.tight);
    const AuAuMLPRuntimeConfig leaderAuAuMLPRuntime = resolveAuAuMLPRuntimeConfig(leaderAuAuStackRuntime.active ? leaderAuAuStackRuntime.mlpMode : cfg.tight);
    std::string auauRuntimeTightModelFile = leaderAuAuRuntime.modelFile;
    std::vector<std::string> auauRuntimeTightFeatures = leaderAuAuRuntime.features;
    std::vector<std::string> auauRuntimeCentDepModelFiles = leaderAuAuRuntime.centModelFiles;

    if (leaderAuAuRuntime.active)
    {
        tightPhotonNode = "PHOTONCLUSTER_CEMC";
        if (!leaderAuAuRuntime.centEdges.empty())
        {
            for (std::size_t i = 0; i + 1 < leaderAuAuRuntime.centEdges.size(); ++i)
            {
                std::ostringstream scoreName;
                scoreName << "auau_tight_bdt_score_cent_"
                          << leaderAuAuRuntime.centEdges[i] << "_" << leaderAuAuRuntime.centEdges[i + 1];
                auauCentDepScoreNames.push_back(scoreName.str());
            }
        }
        if (leaderAuAuRuntime.usesBuilderScore)
        {
            photonBuilder->add_named_bdt_score("auau_tight_bdt_score",
                                               leaderAuAuRuntime.modelFile,
                                               leaderAuAuRuntime.features,
                                               5.0f,
                                               80.0f,
                                               static_cast<float>(cfg.photon_eta_abs_max));
        }
    }
    if (leaderAuAuMLPRuntime.active)
    {
        tightPhotonNode = "PHOTONCLUSTER_CEMC";
    }
    if (leaderAuAuStackRuntime.active)
    {
        tightPhotonNode = "PHOTONCLUSTER_CEMC";
    }

    if (useSamePhotonBDTScores)
    {
        se->registerSubsystem(photonBuilder);
    }
    
    if (vlevel > 0)
    {
        Dl_info pcbInfo{};
        if (dladdr((void*)&typeid(PhotonClusterBuilder), &pcbInfo) && pcbInfo.dli_fname)
            std::cout << "[DBG] PhotonClusterBuilder RTTI from: " << pcbInfo.dli_fname << "\n";
        else
            std::cout << "[DBG] PhotonClusterBuilder RTTI probe: dladdr failed\n";
        
        std::cout << "[DBG] PhotonClusterBuilder vzCut config: use="
        << (cfg.use_vz_cut ? "true" : "false")
        << " vz_cut_cm=" << photonBuilderVzCutCm
        << " | ppg12RecoVertex="
        << (ppg12PhotonYieldPPSim ? "GlobalVertexMap::MBD" : "default reco")
        << " | ppg12TruthVertexRole=SI/DI weight only"
        << " | cemcShapeDiagnosticVariant=" << cemcShowerShapeDiagnosticVariant
        << " | cemcShapeEnergySource=" << resolvedCEMCShapeEnergySource
        << " | shapeTowerAcceptance=" << resolvedCEMCShapeTowerAcceptance
        << " | rawTowermapCEMCShapes=" << (useRawClusterTowermapForCEMCShapes ? "true" : "false")
        << " | shapeTowerMinEGeV=" << resolvedCEMCShapeTowerMinGeV
        << " | isAuAuLike=" << (isAuAuLike ? "true" : "false")
        << " | isSimEmbedded=" << (isSimEmbedded ? "true" : "false")
        << " | photonBuilderIsAuAu=" << (photonBuilderIsAuAu ? "true" : "false")
        << " | clusterUEpipeline=" << cfg.clusterUEpipeline
        << " | inputClusterNode=" << photonInputClusterNode
        << " | preselectionVariant=" << cfg.preselection
        << " | tightVariant=" << cfg.tight
        << " | preselectionPhotonNode=" << preselectionPhotonNode
        << " | tightPhotonNode=" << tightPhotonNode << "\n";
    }
    
    auto fmtDouble = [](double x) -> std::string
    {
        std::ostringstream os;
        os << std::setprecision(17) << x;
        return os.str();
    };
    auto joinStrings = [](const std::vector<std::string>& vals) -> std::string
    {
        std::ostringstream os;
        for (std::size_t i = 0; i < vals.size(); ++i)
        {
            if (i) os << ',';
            os << vals[i];
        }
        return os.str();
    };
    auto joinInts = [](const std::vector<int>& vals) -> std::string
    {
        std::ostringstream os;
        for (std::size_t i = 0; i < vals.size(); ++i)
        {
            if (i) os << ',';
            os << vals[i];
        }
        return os.str();
    };
    auto joinDoubles = [](const std::vector<double>& vals) -> std::string
    {
        std::ostringstream os;
        for (std::size_t i = 0; i < vals.size(); ++i)
        {
            if (i) os << ',';
            os << vals[i];
        }
        return os.str();
    };
    auto envOrDefault = [](const char* key, const std::string& fallback) -> std::string
    {
        const char* raw = std::getenv(key);
        return raw ? std::string(raw) : fallback;
    };
    
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PRESELECTION_VARIANT",
                                               "RJ_PRESELECTION_VARIANT",
                                               cfg.preselection));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_TIGHT_VARIANT",
                                               "RJ_TIGHT_VARIANT",
                                               cfg.tight));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_NONTIGHT_VARIANT",
                                               "RJ_NONTIGHT_VARIANT",
                                               cfg.nonTight));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PRESELECTION_PHOTON_NODE",
                                               "RJ_PRESELECTION_PHOTON_NODE",
                                               preselectionPhotonNode));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_TIGHT_PHOTON_NODE",
                                               "RJ_TIGHT_PHOTON_NODE",
                                               tightPhotonNode));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_NPB_CUT",
                                               "RJ_NPB_CUT",
                                               fmtDouble(cfg.npb_cut)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NPB_CUT",
                                               "RJ_AUAU_NPB_CUT",
                                               fmtDouble(cfg.auau_npb_cut)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_TIGHT_BDT_MIN_INTERCEPT",
                                               "RJ_TIGHT_BDT_MIN_INTERCEPT",
                                               fmtDouble(cfg.tight_bdt_min_intercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_TIGHT_BDT_MIN_SLOPE",
                                               "RJ_TIGHT_BDT_MIN_SLOPE",
                                               fmtDouble(cfg.tight_bdt_min_slope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_TIGHT_BDT_MAX",
                                               "RJ_TIGHT_BDT_MAX",
                                               fmtDouble(cfg.tight_bdt_max)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_NONTIGHT_BDT_MIN_INTERCEPT",
                                               "RJ_NONTIGHT_BDT_MIN_INTERCEPT",
                                               fmtDouble(cfg.nontight_bdt_min_intercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_NONTIGHT_BDT_MIN_SLOPE",
                                               "RJ_NONTIGHT_BDT_MIN_SLOPE",
                                               fmtDouble(cfg.nontight_bdt_min_slope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_NONTIGHT_BDT_MAX_INTERCEPT",
                                               "RJ_NONTIGHT_BDT_MAX_INTERCEPT",
                                               fmtDouble(cfg.nontight_bdt_max_intercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_NONTIGHT_BDT_MAX_SLOPE",
                                               "RJ_NONTIGHT_BDT_MAX_SLOPE",
                                               fmtDouble(cfg.nontight_bdt_max_slope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_MIN_INTERCEPT",
                                               "RJ_AUAU_TIGHT_BDT_MIN_INTERCEPT",
                                               fmtDouble(cfg.auau_tight_bdt_min_intercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_MIN_SLOPE",
                                               "RJ_AUAU_TIGHT_BDT_MIN_SLOPE",
                                               fmtDouble(cfg.auau_tight_bdt_min_slope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_MAX",
                                               "RJ_AUAU_TIGHT_BDT_MAX",
                                               fmtDouble(cfg.auau_tight_bdt_max)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_APPLY_PT_MIN",
                                               "RJ_AUAU_TIGHT_BDT_APPLY_PT_MIN",
                                               fmtDouble(cfg.auau_tight_bdt_apply_pt_min)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_APPLY_PT_MAX",
                                               "RJ_AUAU_TIGHT_BDT_APPLY_PT_MAX",
                                               fmtDouble(cfg.auau_tight_bdt_apply_pt_max)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_MODEL_FILE",
                                               "RJ_AUAU_TIGHT_BDT_MODEL_FILE",
                                               auauRuntimeTightModelFile));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_SCORING_MODE",
                                               "RJ_AUAU_TIGHT_BDT_SCORING_MODE",
                                               leaderAuAuStackRuntime.active ? leaderAuAuStackRuntime.bdtMode : cfg.tight));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_FEATURES",
                                               "RJ_AUAU_TIGHT_BDT_FEATURES",
                                               joinStrings(auauRuntimeTightFeatures)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_CENTDEP_MODEL_FILES",
                                               "RJ_AUAU_TIGHT_BDT_CENTDEP_MODEL_FILES",
                                               joinStrings(auauRuntimeCentDepModelFiles)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_CENTDEP_EDGES",
                                               "RJ_AUAU_TIGHT_BDT_CENTDEP_EDGES",
                                               joinInts(cfg.centrality_edges)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_CENTDEP_SCORE_NAMES",
                                               "RJ_AUAU_TIGHT_BDT_CENTDEP_SCORE_NAMES",
                                               joinStrings(auauCentDepScoreNames)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_MLP_MODEL_FILE",
                                               "RJ_AUAU_TIGHT_MLP_MODEL_FILE",
                                               leaderAuAuMLPRuntime.modelFile));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_MLP_SCORING_MODE",
                                               "RJ_AUAU_TIGHT_MLP_SCORING_MODE",
                                               leaderAuAuStackRuntime.active ? leaderAuAuStackRuntime.mlpMode : cfg.tight));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_MLP_MIN_INTERCEPT",
                                               "RJ_AUAU_TIGHT_MLP_MIN_INTERCEPT",
                                               fmtDouble(cfg.auau_tight_mlp_min_intercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_MLP_MIN_SLOPE",
                                               "RJ_AUAU_TIGHT_MLP_MIN_SLOPE",
                                               fmtDouble(cfg.auau_tight_mlp_min_slope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_MLP_MAX",
                                               "RJ_AUAU_TIGHT_MLP_MAX",
                                               fmtDouble(cfg.auau_tight_mlp_max)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_MLP_APPLY_PT_MIN",
                                               "RJ_AUAU_TIGHT_MLP_APPLY_PT_MIN",
                                               fmtDouble(cfg.auau_tight_mlp_apply_pt_min)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_MLP_APPLY_PT_MAX",
                                               "RJ_AUAU_TIGHT_MLP_APPLY_PT_MAX",
                                               fmtDouble(cfg.auau_tight_mlp_apply_pt_max)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_MLP_WORKING_POINT_ENTRIES",
                                               "RJ_AUAU_TIGHT_MLP_WORKING_POINT_ENTRIES",
                                               joinStrings(cfg.auau_tight_mlp_working_point_entries)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_MLP_MIN_INTERCEPT",
                                               "RJ_AUAU_NONTIGHT_MLP_MIN_INTERCEPT",
                                               fmtDouble(cfg.auau_nontight_mlp_min_intercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_MLP_MIN_SLOPE",
                                               "RJ_AUAU_NONTIGHT_MLP_MIN_SLOPE",
                                               fmtDouble(cfg.auau_nontight_mlp_min_slope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_MLP_MAX_INTERCEPT",
                                               "RJ_AUAU_NONTIGHT_MLP_MAX_INTERCEPT",
                                               fmtDouble(cfg.auau_nontight_mlp_max_intercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_MLP_MAX_SLOPE",
                                               "RJ_AUAU_NONTIGHT_MLP_MAX_SLOPE",
                                               fmtDouble(cfg.auau_nontight_mlp_max_slope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_MLP_STACK_MODEL_FILE",
                                               "RJ_AUAU_TIGHT_BDT_MLP_STACK_MODEL_FILE",
                                               leaderAuAuStackRuntime.modelFile));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_MLP_STACK_MIN_INTERCEPT",
                                               "RJ_AUAU_TIGHT_BDT_MLP_STACK_MIN_INTERCEPT",
                                               fmtDouble(leaderAuAuStackRuntime.minIntercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_MLP_STACK_MIN_SLOPE",
                                               "RJ_AUAU_TIGHT_BDT_MLP_STACK_MIN_SLOPE",
                                               fmtDouble(leaderAuAuStackRuntime.minSlope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_MLP_STACK_MAX",
                                               "RJ_AUAU_TIGHT_BDT_MLP_STACK_MAX",
                                               fmtDouble(leaderAuAuStackRuntime.maxScore)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_MLP_STACK_APPLY_PT_MIN",
                                               "RJ_AUAU_TIGHT_BDT_MLP_STACK_APPLY_PT_MIN",
                                               fmtDouble(leaderAuAuStackRuntime.applyPtMin)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_MLP_STACK_APPLY_PT_MAX",
                                               "RJ_AUAU_TIGHT_BDT_MLP_STACK_APPLY_PT_MAX",
                                               fmtDouble(leaderAuAuStackRuntime.applyPtMax)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_BDT_MLP_STACK_WORKING_POINT_ENTRIES",
                                               "RJ_AUAU_TIGHT_BDT_MLP_STACK_WORKING_POINT_ENTRIES",
                                               joinStrings(leaderAuAuStackRuntime.workingPointEntries)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_BDT_MLP_STACK_MIN_INTERCEPT",
                                               "RJ_AUAU_NONTIGHT_BDT_MLP_STACK_MIN_INTERCEPT",
                                               fmtDouble(leaderAuAuStackRuntime.nonTightMinIntercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_BDT_MLP_STACK_MIN_SLOPE",
                                               "RJ_AUAU_NONTIGHT_BDT_MLP_STACK_MIN_SLOPE",
                                               fmtDouble(leaderAuAuStackRuntime.nonTightMinSlope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_BDT_MLP_STACK_MAX_INTERCEPT",
                                               "RJ_AUAU_NONTIGHT_BDT_MLP_STACK_MAX_INTERCEPT",
                                               fmtDouble(leaderAuAuStackRuntime.nonTightMaxIntercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_BDT_MLP_STACK_MAX_SLOPE",
                                               "RJ_AUAU_NONTIGHT_BDT_MLP_STACK_MAX_SLOPE",
                                               fmtDouble(leaderAuAuStackRuntime.nonTightMaxSlope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_LOGREG_MODEL_FILE",
                                               "RJ_AUAU_TIGHT_LOGREG_MODEL_FILE",
                                               leaderAuAuLogRegRuntime.modelFile));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_LOGREG_MIN_INTERCEPT",
                                               "RJ_AUAU_TIGHT_LOGREG_MIN_INTERCEPT",
                                               fmtDouble(leaderAuAuLogRegRuntime.minIntercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_LOGREG_MIN_SLOPE",
                                               "RJ_AUAU_TIGHT_LOGREG_MIN_SLOPE",
                                               fmtDouble(leaderAuAuLogRegRuntime.minSlope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_LOGREG_MAX",
                                               "RJ_AUAU_TIGHT_LOGREG_MAX",
                                               fmtDouble(leaderAuAuLogRegRuntime.maxScore)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_LOGREG_APPLY_PT_MIN",
                                               "RJ_AUAU_TIGHT_LOGREG_APPLY_PT_MIN",
                                               fmtDouble(leaderAuAuLogRegRuntime.applyPtMin)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_LOGREG_APPLY_PT_MAX",
                                               "RJ_AUAU_TIGHT_LOGREG_APPLY_PT_MAX",
                                               fmtDouble(leaderAuAuLogRegRuntime.applyPtMax)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_TIGHT_LOGREG_WORKING_POINT_ENTRIES",
                                               "RJ_AUAU_TIGHT_LOGREG_WORKING_POINT_ENTRIES",
                                               joinStrings(leaderAuAuLogRegRuntime.workingPointEntries)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_LOGREG_MIN_INTERCEPT",
                                               "RJ_AUAU_NONTIGHT_LOGREG_MIN_INTERCEPT",
                                               fmtDouble(leaderAuAuLogRegRuntime.nonTightMinIntercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_LOGREG_MIN_SLOPE",
                                               "RJ_AUAU_NONTIGHT_LOGREG_MIN_SLOPE",
                                               fmtDouble(leaderAuAuLogRegRuntime.nonTightMinSlope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_LOGREG_MAX_INTERCEPT",
                                               "RJ_AUAU_NONTIGHT_LOGREG_MAX_INTERCEPT",
                                               fmtDouble(leaderAuAuLogRegRuntime.nonTightMaxIntercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_LOGREG_MAX_SLOPE",
                                               "RJ_AUAU_NONTIGHT_LOGREG_MAX_SLOPE",
                                               fmtDouble(leaderAuAuLogRegRuntime.nonTightMaxSlope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_BDT_MIN_INTERCEPT",
                                               "RJ_AUAU_NONTIGHT_BDT_MIN_INTERCEPT",
                                               fmtDouble(cfg.auau_nontight_bdt_min_intercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_BDT_MIN_SLOPE",
                                               "RJ_AUAU_NONTIGHT_BDT_MIN_SLOPE",
                                               fmtDouble(cfg.auau_nontight_bdt_min_slope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_BDT_MAX_INTERCEPT",
                                               "RJ_AUAU_NONTIGHT_BDT_MAX_INTERCEPT",
                                               fmtDouble(cfg.auau_nontight_bdt_max_intercept)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_BDT_MAX_SLOPE",
                                               "RJ_AUAU_NONTIGHT_BDT_MAX_SLOPE",
                                               fmtDouble(cfg.auau_nontight_bdt_max_slope)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_BDT_SIDEBAND_MODE",
                                               "RJ_AUAU_NONTIGHT_BDT_SIDEBAND_MODE",
                                               cfg.auau_nontight_bdt_sideband_mode));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_BDT_RELATIVE_MIN_OFFSET",
                                               "RJ_AUAU_NONTIGHT_BDT_RELATIVE_MIN_OFFSET",
                                               fmtDouble(cfg.auau_nontight_bdt_relative_min_offset)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NONTIGHT_BDT_RELATIVE_MAX_OFFSET",
                                               "RJ_AUAU_NONTIGHT_BDT_RELATIVE_MAX_OFFSET",
                                               fmtDouble(cfg.auau_nontight_bdt_relative_max_offset)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_BDT_TRAINING_TREE",
                                               "RJ_AUAU_BDT_TRAINING_TREE",
                                               cfg.auau_bdt_training_tree ? "true" : "false"));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_BDT_TRAINING_TREE_MAX_ENTRIES",
                                               "RJ_AUAU_BDT_TRAINING_TREE_MAX_ENTRIES",
                                               std::to_string(cfg.auau_bdt_training_tree_max_entries)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_BDT_NPB_DATA_TAGGING",
                                               "RJ_AUAU_BDT_NPB_DATA_TAGGING",
                                               cfg.auau_bdt_npb_data_tagging ? "true" : "false"));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PP_PHOTONID_EXTRACT_ONLY",
                                               "RJ_PP_PHOTONID_EXTRACT_ONLY",
                                               envOrDefault("RJ_PP_PHOTONID_EXTRACT_ONLY",
                                                            cfg.pp_photonid_extract_only ? "true" : "false")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PP_PHOTONID_TRAINING_TREE",
                                               "RJ_PP_PHOTONID_TRAINING_TREE",
                                               envOrDefault("RJ_PP_PHOTONID_TRAINING_TREE",
                                                            cfg.pp_photonid_training_tree ? "true" : "false")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES",
                                               "RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES",
                                               envOrDefault("RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES",
                                                            std::to_string(cfg.pp_photonid_training_tree_max_entries))));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PP_PHOTONID_SOURCE_ROLE",
                                               "RJ_PP_PHOTONID_SOURCE_ROLE",
                                               envOrDefault("RJ_PP_PHOTONID_SOURCE_ROLE",
                                                            cfg.pp_photonid_source_role)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PP_PHOTONID_PPG12_FILTER",
                                               "RJ_PP_PHOTONID_PPG12_FILTER",
                                               envOrDefault("RJ_PP_PHOTONID_PPG12_FILTER",
                                                            cfg.pp_photonid_ppg12_filter ? "true" : "false")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PP_PHOTONID_REQUIRE_PRESELECTION",
                                               "RJ_PP_PHOTONID_REQUIRE_PRESELECTION",
                                               envOrDefault("RJ_PP_PHOTONID_REQUIRE_PRESELECTION",
                                                            cfg.pp_photonid_require_preselection ? "true" : "false")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_TABLE_QA",
                                               "RJ_PPG12_TABLE_QA",
                                               envOrDefault("RJ_PPG12_TABLE_QA", "0")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_TABLE_QA_MC_ISO_SCALE",
                                               "RJ_PPG12_TABLE_QA_MC_ISO_SCALE",
                                               envOrDefault("RJ_PPG12_TABLE_QA_MC_ISO_SCALE", "1.2")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_TABLE_QA_MC_ISO_SHIFT",
                                               "RJ_PPG12_TABLE_QA_MC_ISO_SHIFT",
                                               envOrDefault("RJ_PPG12_TABLE_QA_MC_ISO_SHIFT", "0.2")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_TABLE_QA_NPB_DATA_TAGGING",
                                               "RJ_PPG12_TABLE_QA_NPB_DATA_TAGGING",
                                               envOrDefault("RJ_PPG12_TABLE_QA_NPB_DATA_TAGGING", "1")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_TABLE_QA_NPB_TIME_SAMPLE_NS",
                                               "RJ_PPG12_TABLE_QA_NPB_TIME_SAMPLE_NS",
                                               envOrDefault("RJ_PPG12_TABLE_QA_NPB_TIME_SAMPLE_NS", "17.6")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_TABLE_QA_NPB_DELTA_T_CUT",
                                               "RJ_PPG12_TABLE_QA_NPB_DELTA_T_CUT",
                                               envOrDefault("RJ_PPG12_TABLE_QA_NPB_DELTA_T_CUT", "-5.0")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_TABLE_QA_NPB_WETA_MIN",
                                               "RJ_PPG12_TABLE_QA_NPB_WETA_MIN",
                                               envOrDefault("RJ_PPG12_TABLE_QA_NPB_WETA_MIN", "0.4")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_TABLE_QA_NPB_AWAY_JET_PT_MIN",
                                               "RJ_PPG12_TABLE_QA_NPB_AWAY_JET_PT_MIN",
                                               envOrDefault("RJ_PPG12_TABLE_QA_NPB_AWAY_JET_PT_MIN", "5.0")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_TABLE_QA_NPB_AWAY_JET_DPHI_MIN",
                                               "RJ_PPG12_TABLE_QA_NPB_AWAY_JET_DPHI_MIN",
                                               envOrDefault("RJ_PPG12_TABLE_QA_NPB_AWAY_JET_DPHI_MIN", "1.5707963267948966")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_TABLE_QA_MBD_T0_CORRECTION_FILE",
                                               "RJ_PPG12_TABLE_QA_MBD_T0_CORRECTION_FILE",
                                               envOrDefault("RJ_PPG12_TABLE_QA_MBD_T0_CORRECTION_FILE", "")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_FIG11_SB_DIAGNOSTIC",
                                               "RJ_PPG12_FIG11_SB_DIAGNOSTIC",
                                               envOrDefault("RJ_PPG12_FIG11_SB_DIAGNOSTIC", "0")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_FIG13_PARITY_QA",
                                               "RJ_PPG12_FIG13_PARITY_QA",
                                               envOrDefault("RJ_PPG12_FIG13_PARITY_QA", "0")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_FIG7_TRIGGER_DIAGNOSTIC",
                                               "RJ_PPG12_FIG7_TRIGGER_DIAGNOSTIC",
                                               envOrDefault("RJ_PPG12_FIG7_TRIGGER_DIAGNOSTIC", "0")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_PPG12_FIG13_BIT30_DIAGNOSTIC",
                                               "RJ_PPG12_FIG13_BIT30_DIAGNOSTIC",
                                               envOrDefault("RJ_PPG12_FIG13_BIT30_DIAGNOSTIC", "0")));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NPB_TAG_DELTA_T_CUT",
                                               "RJ_AUAU_NPB_TAG_DELTA_T_CUT",
                                               fmtDouble(cfg.auau_npb_tag_delta_t_cut)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NPB_TAG_WETA_MIN",
                                               "RJ_AUAU_NPB_TAG_WETA_MIN",
                                               fmtDouble(cfg.auau_npb_tag_weta_min)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NPB_TAG_AWAY_JET_PT_MIN",
                                               "RJ_AUAU_NPB_TAG_AWAY_JET_PT_MIN",
                                               fmtDouble(cfg.auau_npb_tag_away_jet_pt_min)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NPB_TAG_AWAY_JET_DPHI_MIN",
                                               "RJ_AUAU_NPB_TAG_AWAY_JET_DPHI_MIN",
                                               fmtDouble(cfg.auau_npb_tag_away_jet_dphi_min)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NPB_TAG_TIME_SAMPLE_NS",
                                               "RJ_AUAU_NPB_TAG_TIME_SAMPLE_NS",
                                               fmtDouble(cfg.auau_npb_tag_time_sample_ns)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_AUAU_NPB_MBD_T0_OFFSET",
                                               "RJ_AUAU_NPB_MBD_T0_OFFSET",
                                               fmtDouble(cfg.auau_npb_mbd_t0_offset)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_JET_ML_TRAINING_TREE",
                                               "RJ_JET_ML_TRAINING_TREE",
                                               cfg.jet_ml_training_tree ? "true" : "false"));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_JET_ML_TRAINING_TREE_MAX_ENTRIES",
                                               "RJ_JET_ML_TRAINING_TREE_MAX_ENTRIES",
                                               std::to_string(cfg.jet_ml_training_tree_max_entries)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_JET_ML_CORRECTION_ENABLED",
                                               "RJ_JET_ML_CORRECTION_ENABLED",
                                               cfg.jet_ml_correction_enabled ? "true" : "false"));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_JET_ML_MODEL_FILE",
                                               "RJ_JET_ML_MODEL_FILE",
                                               cfg.jet_ml_model_file));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_JET_ML_FEATURES",
                                               "RJ_JET_ML_FEATURES",
                                               joinStrings(cfg.jet_ml_features)));
    se->registerSubsystem(new ProcessEnvSetter("Env_RJ_JET_ML_USE_AS_NOMINAL",
                                               "RJ_JET_ML_USE_AS_NOMINAL",
                                               cfg.jet_ml_use_as_nominal ? "true" : "false"));
    
    // ------------------------------------------------------------------
    // Apply YAML-driven knobs. In fanout mode each RecoilJets instance gets
    // its own output file and photon-ID triplet, but shares the expensive
    // upstream DST, cluster, photon-score, and jet reconstruction pass.
    // ------------------------------------------------------------------
    auto configureRecoilJetsInstance =
    [&](RecoilJets* recoilJets, const idfanout::Entry& idEntry)
    {
	    recoilJets->setPhotonIDVariants(idEntry.preselection,
	                                    idEntry.tight,
	                                    idEntry.nonTight,
	                                    preselectionPhotonNode,
	                                    tightPhotonNode);
        const AuAuBDTMLPStackRuntimeConfig entryAuAuStackRuntime = resolveAuAuBDTMLPStackRuntimeConfig(idEntry.tight);
        const AuAuBDTRuntimeConfig entryAuAuRuntime = resolveAuAuBDTRuntimeConfig(entryAuAuStackRuntime.active ? entryAuAuStackRuntime.bdtMode : idEntry.tight);
        if (entryAuAuRuntime.active)
        {
            recoilJets->setAuAuTightBDTRuntimeConfig(entryAuAuRuntime.modelFile,
                                                     entryAuAuRuntime.features,
                                                     entryAuAuRuntime.centEdges,
                                                     entryAuAuRuntime.centModelFiles,
                                                     entryAuAuRuntime.ptEdges,
                                                     entryAuAuRuntime.ptModelFiles,
                                                     entryAuAuRuntime.ptCentModelFiles,
                                                     entryAuAuRuntime.ptFallbackModelFile,
                                                     entryAuAuRuntime.ptFallbackCentModelFiles,
                                                     entryAuAuRuntime.ptFallbackMin,
                                                     entryAuAuRuntime.ptFallbackMax,
                                                     entryAuAuRuntime.applyPtMin,
                                                     entryAuAuRuntime.applyPtMax,
                                                     entryAuAuRuntime.workingPointEntries,
                                                     entryAuAuStackRuntime.active ? entryAuAuStackRuntime.bdtMode : idEntry.tight);
        }
        const AuAuMLPRuntimeConfig entryAuAuMLPRuntime = resolveAuAuMLPRuntimeConfig(entryAuAuStackRuntime.active ? entryAuAuStackRuntime.mlpMode : idEntry.tight);
        if (entryAuAuMLPRuntime.active)
        {
            recoilJets->setAuAuTightMLPRuntimeConfig(entryAuAuMLPRuntime.modelFile,
                                                     entryAuAuMLPRuntime.minIntercept,
                                                     entryAuAuMLPRuntime.minSlope,
                                                     entryAuAuMLPRuntime.maxScore,
                                                     entryAuAuMLPRuntime.nonTightMinIntercept,
                                                     entryAuAuMLPRuntime.nonTightMinSlope,
                                                     entryAuAuMLPRuntime.nonTightMaxIntercept,
                                                     entryAuAuMLPRuntime.nonTightMaxSlope,
                                                     entryAuAuMLPRuntime.applyPtMin,
                                                     entryAuAuMLPRuntime.applyPtMax,
                                                     entryAuAuMLPRuntime.workingPointEntries,
                                                     entryAuAuStackRuntime.active ? entryAuAuStackRuntime.mlpMode : idEntry.tight);
        }
        if (entryAuAuStackRuntime.active)
        {
            recoilJets->setAuAuTightBDTMLPStackRuntimeConfig(entryAuAuStackRuntime.modelFile,
                                                            entryAuAuStackRuntime.minIntercept,
                                                            entryAuAuStackRuntime.minSlope,
                                                            entryAuAuStackRuntime.maxScore,
                                                            entryAuAuStackRuntime.nonTightMinIntercept,
                                                            entryAuAuStackRuntime.nonTightMinSlope,
                                                            entryAuAuStackRuntime.nonTightMaxIntercept,
                                                            entryAuAuStackRuntime.nonTightMaxSlope,
                                                            entryAuAuStackRuntime.applyPtMin,
                                                            entryAuAuStackRuntime.applyPtMax,
                                                            entryAuAuStackRuntime.workingPointEntries);
        }
        const AuAuLogRegRuntimeConfig entryAuAuLogRegRuntime = resolveAuAuLogRegRuntimeConfig(idEntry.tight);
        if (entryAuAuLogRegRuntime.active)
        {
            recoilJets->setAuAuTightLogRegRuntimeConfig(entryAuAuLogRegRuntime.modelFile,
                                                        entryAuAuLogRegRuntime.minIntercept,
                                                        entryAuAuLogRegRuntime.minSlope,
                                                        entryAuAuLogRegRuntime.maxScore,
                                                        entryAuAuLogRegRuntime.nonTightMinIntercept,
                                                        entryAuAuLogRegRuntime.nonTightMinSlope,
                                                        entryAuAuLogRegRuntime.nonTightMaxIntercept,
                                                        entryAuAuLogRegRuntime.nonTightMaxSlope,
                                                        entryAuAuLogRegRuntime.applyPtMin,
                                                        entryAuAuLogRegRuntime.applyPtMax,
                                                        entryAuAuLogRegRuntime.workingPointEntries);
        }
	    recoilJets->setPhotonEtaAbsMax(cfg.photon_eta_abs_max);
    recoilJets->setMinJetPt(cfg.jet_pt_min);
    recoilJets->setMinBackToBack(cfg.back_to_back_dphi_min_pi_fraction * M_PI);
    
    recoilJets->setUseVzCut(cfg.use_vz_cut, cfg.vz_cut_cm);
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
    recoilJets->setMinBiasClassifier(cfg.setMinBiasClassifer);
    recoilJets->setCentEdges(cfg.centrality_edges);
    recoilJets->setRequireTowerInfoTruthMatching(cfg.require_towerinfo_truth_matching);
    recoilJets->setVertexReweighting(cfg.vertex_reweight_on_auau,
                                     cfg.vertex_reweight_file_auau,
                                     cfg.vertex_reweight_hist_auau);
    recoilJets->setCentralityReweighting(cfg.centrality_reweight_on,
                                         cfg.centrality_reweight_file,
                                         cfg.centrality_reweight_hist);
#else
    bool ppVertexReweightOn = cfg.vertex_reweight_on_pp;
    std::string ppVertexReweightFile = cfg.vertex_reweight_file_pp;
    std::string ppVertexReweightHist = cfg.vertex_reweight_hist_pp;
    if (const char* env = std::getenv("RJ_PP_VERTEX_REWEIGHT_FILE"))
    {
        const std::string value = detail::trim(std::string(env));
        if (!value.empty())
        {
            ppVertexReweightOn = true;
            ppVertexReweightFile = value;
        }
    }
    if (const char* env = std::getenv("RJ_PP_VERTEX_REWEIGHT_HIST"))
    {
        const std::string value = detail::trim(std::string(env));
        if (!value.empty()) ppVertexReweightHist = value;
    }
    recoilJets->setVertexReweighting(ppVertexReweightOn,
                                     ppVertexReweightFile,
                                     ppVertexReweightHist);
#endif
    recoilJets->setActiveJetRKeys(activeJetRKeys);
    const double entryConeR = idEntry.hasIsoOverride ? idEntry.coneR : cfg.isoConeR;
    const bool entryIsSlidingIso = idEntry.hasIsoOverride ? idEntry.isSlidingIso : cfg.isSlidingIso;
    const double entryFixedGeV = idEntry.hasIsoOverride ? idEntry.fixedGeV : cfg.isoFixed;
    recoilJets->setIsolationWP(cfg.isoA, cfg.isoB, cfg.isoGap, entryConeR, recoilJetsIsoTowerMin, entryFixedGeV);
    recoilJets->setIsSlidingIso(entryIsSlidingIso);
    recoilJets->setTruthIsoMaxGeV(cfg.truthIsoGeV);
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
    if (!cfg.auauCentIsoWP.empty())
    {
        std::vector<RecoilJets::CentIsoWP> wps;
        wps.reserve(cfg.auauCentIsoWP.size());
        for (const auto& e : cfg.auauCentIsoWP)
            wps.push_back({e.aGeV, e.bPerGeV, e.sideGapGeV});
        recoilJets->setCentIsoWPs(wps);
    }
    if (!cfg.auauCentIsoWPR30.empty())
    {
        std::vector<RecoilJets::CentIsoWP> wps;
        wps.reserve(cfg.auauCentIsoWPR30.size());
        for (const auto& e : cfg.auauCentIsoWPR30)
            wps.push_back({e.aGeV, e.bPerGeV, e.sideGapGeV});
        recoilJets->setCentIsoWPsForCone(0.30, wps);
    }
    if (!cfg.auauCentIsoWPR40.empty())
    {
        std::vector<RecoilJets::CentIsoWP> wps;
        wps.reserve(cfg.auauCentIsoWPR40.size());
        for (const auto& e : cfg.auauCentIsoWPR40)
            wps.push_back({e.aGeV, e.bPerGeV, e.sideGapGeV});
        recoilJets->setCentIsoWPsForCone(0.40, wps);
    }
#else
    if (cfg.ppIsoWPR30.configured)
        recoilJets->setPPIsoWPForCone(0.30, cfg.ppIsoWPR30.aGeV, cfg.ppIsoWPR30.bPerGeV, cfg.ppIsoWPR30.sideGapGeV);
    if (cfg.ppIsoWPR40.configured)
        recoilJets->setPPIsoWPForCone(0.40, cfg.ppIsoWPR40.aGeV, cfg.ppIsoWPR40.bPerGeV, cfg.ppIsoWPR40.sideGapGeV);
#endif
    
    recoilJets->setPhotonIDCuts(cfg.pre_e11e33_max,
                                cfg.pre_et1_min,
                                cfg.pre_et1_max,
                                cfg.pre_e32e35_min,
                                cfg.pre_e32e35_max,
                                cfg.pre_weta_max,
                                cfg.tight_w_lo,
                                cfg.tight_w_hi_intercept,
                                cfg.tight_w_hi_slope,
                                cfg.tight_e11e33_min,
                                cfg.tight_e11e33_max,
                                cfg.tight_et1_min,
                                cfg.tight_et1_max,
                                cfg.tight_e32e35_min,
                                cfg.tight_e32e35_max);
    
    recoilJets->setGammaPtBins(cfg.jes3_photon_pt_bins);
    recoilJets->setPhoMatchDRMax(cfg.pho_dr_max);
    recoilJets->setJetMatchDRMax(cfg.jet_dr_max);
    
    recoilJets->setUnfoldRecoPhotonPtBins(cfg.unfold_reco_photon_pt_bins);
    recoilJets->setUnfoldTruthPhotonPtBins(cfg.unfold_truth_photon_pt_bins);
    recoilJets->setUnfoldJetPtBins(unfoldJetPtEdges);
    recoilJets->setUnfoldXJBins(cfg.unfold_xj_bins);
    recoilJets->setLeadingResponseFamilyLabel(cfg.leading_response_family);
    
    recoilJets->enablePi0Analysis(cfg.doPi0Analysis);
    std::string stampedYaml = idfanout::YAMLForEntry(cfg.yamlText, idEntry);
    idfanout::ReplaceOrAppendScalar(stampedYaml,
                                    "cemc_shower_shape_energy_source",
                                    resolvedCEMCShapeEnergySource);
    idfanout::ReplaceOrAppendScalar(stampedYaml,
                                    "cemc_shower_shape_raw_cluster_towermap",
                                    useRawClusterTowermapForCEMCShapes ? "true" : "false");
    idfanout::ReplaceOrAppendScalar(stampedYaml,
                                    "cemc_shower_shape_tower_acceptance",
                                    resolvedCEMCShapeTowerAcceptance);
    idfanout::ReplaceOrAppendScalar(stampedYaml,
                                    "cemc_shower_shape_tower_min_energy_gev",
                                    detail::fmt(resolvedCEMCShapeTowerMinGeV, 3));
    idfanout::ReplaceOrAppendScalar(stampedYaml,
                                    "cemc_shower_shape_diagnostic_variant",
                                    cemcShowerShapeDiagnosticVariant);
    idfanout::ReplaceOrAppendScalar(stampedYaml,
                                    "ppg12_di_archived_reco_chain",
                                    usePPG12ArchivedDIG4OnlyReco ? "true" : "false");
    if (usePPG12ArchivedDIG4OnlyReco)
    {
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_di_reconstruction_graph",
            "deployed_run28_g4hits_truthjets_current_release_abi");
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_di_sim_sample", env_lower("RJ_SIM_SAMPLE"));
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_di_truth_jets_mode",
            env_lower("RJ_TRUTH_JETS_MODE", "auto"));
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_di_pedestal_file",
            ppg12PedestalFileForProvenance);
    }
    if (ppg12ClosureCanary)
    {
        idfanout::ReplaceOrAppendScalar(stampedYaml,
                                        "ppg12_closure_canary", "true");
        idfanout::ReplaceOrAppendScalar(stampedYaml,
                                        "ppg12_closure_canary_id",
                                        ppg12ClosureCanaryId);
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_closure_rng_contract",
            "historical_phrandomseed_fifo_replay_v2");
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_closure_randomseed", "absent");
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_closure_seed_replay_sequence",
            ppg12ReplaySeedSequence);
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_closure_pedestal_sequence",
            std::to_string(ppg12NaturalPedestalSequence));
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_closure_pedestal_expected_sequence",
            std::to_string(ppg12ExpectedPedestalSequence));
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_closure_pedestal_file",
            ppg12PedestalFileForProvenance);
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_closure_pedestal_natural_sequence",
            std::to_string(ppg12NaturalPedestalSequence));
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_closure_phrandomseed_call_consumed", "true");
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "ppg12_closure_rebuild_input_mode",
            usePPG12PPSimG4OnlyInput ? "g4_only" : "four_lane");
    }
    if (ppg12DINeutralityCanary)
    {
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "replay_foundation_di_neutrality_canary", "true");
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "replay_foundation_di_neutrality_canary_id",
            ppg12DINeutralityCanaryId);
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "replay_foundation_di_neutrality_rng_contract",
            "historical_phrandomseed_fifo_replay_v2");
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "replay_foundation_di_neutrality_seed_sequence",
            ppg12ReplaySeedSequence);
        idfanout::ReplaceOrAppendScalar(
            stampedYaml, "replay_foundation_di_neutrality_pedestal_sequence",
            std::to_string(ppg12NaturalPedestalSequence));
    }
    if (const char* jetPtRaw = std::getenv("RJ_INTERNAL_JET_PT_MINS"))
    {
        const char* disableRaw = std::getenv("RJ_DISABLE_JET_PT_INTERNALIZATION");
        const bool disabled = disableRaw &&
                              (std::string(disableRaw) == "1" ||
                               std::string(disableRaw) == "true" ||
                               std::string(disableRaw) == "TRUE" ||
                               std::string(disableRaw) == "yes" ||
                               std::string(disableRaw) == "YES" ||
                               std::string(disableRaw) == "on" ||
                               std::string(disableRaw) == "ON");
        if (!disabled && *jetPtRaw)
        {
            std::string text(jetPtRaw);
            for (char& c : text)
            {
                if (c == ',' || c == ';' || c == ':') c = ' ';
            }
            std::stringstream ss(text);
            std::string token;
            std::ostringstream list;
            bool first = true;
            list << "[";
            while (ss >> token)
            {
                if (!first) list << ", ";
                first = false;
                list << token;
            }
            list << "]";
            if (!first)
            {
                idfanout::ReplaceOrAppendScalar(stampedYaml,
                                                "internal_jet_pt_min",
                                                list.str());
            }
        }
    }
    if (const char* dphiRaw = std::getenv("RJ_INTERNAL_DPHI_PI_FRACTIONS"))
    {
        const char* disableRaw = std::getenv("RJ_DISABLE_DPHI_INTERNALIZATION");
        const bool disabled = disableRaw &&
                              (std::string(disableRaw) == "1" ||
                               std::string(disableRaw) == "true" ||
                               std::string(disableRaw) == "TRUE" ||
                               std::string(disableRaw) == "yes" ||
                               std::string(disableRaw) == "YES" ||
                               std::string(disableRaw) == "on" ||
                               std::string(disableRaw) == "ON");
        if (!disabled && *dphiRaw)
        {
            std::string text(dphiRaw);
            for (char& c : text)
            {
                if (c == ',' || c == ';' || c == ':') c = ' ';
            }
            std::stringstream ss(text);
            std::string token;
            std::ostringstream list;
            bool first = true;
            list << "[";
            while (ss >> token)
            {
                if (!first) list << ", ";
                first = false;
                list << token;
            }
            list << "]";
            if (!first)
            {
                idfanout::ReplaceOrAppendScalar(stampedYaml,
                                                "internal_back_to_back_dphi_min_pi_fraction",
                                                list.str());
            }
        }
    }
    if (const char* isoViewsRaw = std::getenv("RJ_INTERNAL_ISO_VIEWS"))
    {
        const char* disableRaw = std::getenv("RJ_DISABLE_ISO_CONE_INTERNALIZATION");
        const bool disabled = disableRaw &&
                              (std::string(disableRaw) == "1" ||
                               std::string(disableRaw) == "true" ||
                               std::string(disableRaw) == "TRUE" ||
                               std::string(disableRaw) == "yes" ||
                               std::string(disableRaw) == "YES" ||
                               std::string(disableRaw) == "on" ||
                               std::string(disableRaw) == "ON");
        if (!disabled && *isoViewsRaw)
        {
            idfanout::ReplaceOrAppendScalar(stampedYaml,
                                            "internal_iso_cone_views",
                                            std::string(isoViewsRaw));
        }
    }
    recoilJets->setAnalysisConfigYAML(stampedYaml, "analysis_config.yaml");
    
    if (vlevel > 0)
    {
        std::cout << "[CFG] Applied to RecoilJets:"
        << " etaAbsMax=" << cfg.photon_eta_abs_max
        << " jetPtMin=" << cfg.jet_pt_min
        << " dphiMin(rad)=" << (cfg.back_to_back_dphi_min_pi_fraction * M_PI)
        << " useVzCut=" << (cfg.use_vz_cut ? "true" : "false")
        << " vzCut=" << cfg.vz_cut_cm
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
        << " minBiasClassifierGate=" << (cfg.setMinBiasClassifer ? "true" : "false")
        << " centEdges_n=" << cfg.centrality_edges.size()
#endif
        << " phoDR=" << cfg.pho_dr_max
        << " jetDR=" << cfg.jet_dr_max
        << " isoTowerMinApplied=" << recoilJetsIsoTowerMin
        << " coneR=" << entryConeR
        << " fixedGeV=" << entryFixedGeV
        << " isSlidingIso=" << (entryIsSlidingIso ? "true" : "false")
        << (isAuAuLike ? " (AuAu forced no isolation tower floor)" : "")
        << " preselection=" << idEntry.preselection
        << " tight=" << idEntry.tight
        << " nonTight=" << idEntry.nonTight
        << " output=" << idEntry.outRoot
        << "\n";
        std::cout << "[CFG] photon nodes:"
        << " base=PHOTONCLUSTER_CEMC"
        << " preselectionNode=" << preselectionPhotonNode
        << " tightNode=" << tightPhotonNode << "\n";
        
        if (yamlcfg::PreselectionUsesNPB(idEntry.preselection))
        {
            std::cout << "[CFG] NPB preselection:"
            << " variant=" << idEntry.preselection
            << " model=" << cfg.npb_model_file
            << " cut=(score > " << cfg.npb_cut << ")"
            << " features_n=" << cfg.npb_features.size() << "\n";
        }

        if (yamlcfg::PreselectionUsesAuAuNPB(idEntry.preselection))
        {
            std::cout << "[CFG] AuAu NPB preselection:"
            << " variant=" << idEntry.preselection
            << " model=" << cfg.auau_npb_model_file
            << " cut=(score > " << cfg.auau_npb_cut << ")"
            << " features_n=" << cfg.auau_npb_features.size() << "\n";
        }
        
        if (idEntry.tight == "newPPG12")
        {
            std::cout << "[CFG] tight BDT:"
            << " model=" << cfg.tight_bdt_model_file
            << " cut=(score > " << cfg.tight_bdt_min_slope << " * ET + " << cfg.tight_bdt_min_intercept
            << ")"
            << " features_n=" << cfg.tight_bdt_features.size();
            if (idEntry.nonTight == "newPPG12")
            {
                std::cout << " nonTight=(" << cfg.nontight_bdt_min_slope << " * ET + " << cfg.nontight_bdt_min_intercept
                          << ", " << cfg.nontight_bdt_max_slope << " * ET + " << cfg.nontight_bdt_max_intercept << ")";
            }
            std::cout << "\n";
        }

        if (yamlcfg::IsAuAuTightBDTMode(idEntry.tight))
        {
            const AuAuBDTRuntimeConfig summaryAuAuRuntime = resolveAuAuBDTRuntimeConfig(idEntry.tight);
            std::cout << "[CFG] AuAu tight BDT:"
            << " mode=" << idEntry.tight;
            if (!summaryAuAuRuntime.modelFile.empty())
                std::cout << " model=" << summaryAuAuRuntime.modelFile;
            if (!summaryAuAuRuntime.centModelFiles.empty())
                std::cout << " cent_models=" << summaryAuAuRuntime.centModelFiles.size()
                          << " cent_edges=[" << joinInts(summaryAuAuRuntime.centEdges) << "]";
            if (!summaryAuAuRuntime.ptModelFiles.empty())
                std::cout << " pt_models=" << summaryAuAuRuntime.ptModelFiles.size()
                          << " pt_edges=[" << joinDoubles(summaryAuAuRuntime.ptEdges) << "]"
                          << " fallback=" << summaryAuAuRuntime.ptFallbackModelFile;
            if (!summaryAuAuRuntime.ptCentModelFiles.empty())
                std::cout << " pt_cent_models=" << summaryAuAuRuntime.ptCentModelFiles.size()
                          << " cent_edges=[" << joinInts(summaryAuAuRuntime.centEdges) << "]"
                          << " pt_edges=[" << joinDoubles(summaryAuAuRuntime.ptEdges) << "]"
                          << " fallback_cent_models=" << summaryAuAuRuntime.ptFallbackCentModelFiles.size();
            else
                std::cout << " features_n=" << summaryAuAuRuntime.features.size();
            std::cout
            << " cut=(score > " << cfg.auau_tight_bdt_min_slope << " * ET + " << cfg.auau_tight_bdt_min_intercept
            << " && score < " << cfg.auau_tight_bdt_max << ")";
            if (std::isfinite(summaryAuAuRuntime.applyPtMin) || std::isfinite(summaryAuAuRuntime.applyPtMax))
            {
                std::cout << " apply_pt=["
                          << (std::isfinite(summaryAuAuRuntime.applyPtMin) ? fmtDouble(summaryAuAuRuntime.applyPtMin) : "-inf")
                          << ","
                          << (std::isfinite(summaryAuAuRuntime.applyPtMax) ? fmtDouble(summaryAuAuRuntime.applyPtMax) : "+inf")
                          << ")";
            }
            if (idEntry.nonTight == "auauBDTSideband")
            {
                std::string sidebandModeLower = detail::trim(cfg.auau_nontight_bdt_sideband_mode);
                std::transform(sidebandModeLower.begin(), sidebandModeLower.end(), sidebandModeLower.begin(),
                               [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
                if (sidebandModeLower == "relativetotight" || sidebandModeLower == "relative_to_tight" ||
                    sidebandModeLower == "tightrelative" || sidebandModeLower == "tight_relative")
                {
                    std::cout << " nonTight=relativeToTight"
                              << " tight+[" << cfg.auau_nontight_bdt_relative_min_offset
                              << "," << cfg.auau_nontight_bdt_relative_max_offset << "]";
                }
                else
                {
                    std::cout << " nonTight=etLinear("
                              << cfg.auau_nontight_bdt_min_slope << " * ET + " << cfg.auau_nontight_bdt_min_intercept
                              << ", " << cfg.auau_nontight_bdt_max_slope << " * ET + " << cfg.auau_nontight_bdt_max_intercept << ")";
                }
            }
            else
            {
                std::cout << " nonTight=complement";
            }
            std::cout << "\n";
        }

        if (yamlcfg::IsAuAuTightMLPMode(idEntry.tight))
        {
            const AuAuMLPRuntimeConfig summaryAuAuMLPRuntime = resolveAuAuMLPRuntimeConfig(idEntry.tight);
            std::cout << "[CFG] AuAu tight MLP:"
                      << " mode=" << idEntry.tight
                      << " model=" << summaryAuAuMLPRuntime.modelFile
                      << " cut=(score > " << cfg.auau_tight_mlp_min_slope << " * ET + " << cfg.auau_tight_mlp_min_intercept
                      << " && score < " << cfg.auau_tight_mlp_max << ")";
            if (std::isfinite(summaryAuAuMLPRuntime.applyPtMin) || std::isfinite(summaryAuAuMLPRuntime.applyPtMax))
            {
                std::cout << " apply_pt=["
                          << (std::isfinite(summaryAuAuMLPRuntime.applyPtMin) ? fmtDouble(summaryAuAuMLPRuntime.applyPtMin) : "-inf")
                          << ","
                          << (std::isfinite(summaryAuAuMLPRuntime.applyPtMax) ? fmtDouble(summaryAuAuMLPRuntime.applyPtMax) : "+inf")
                          << ")";
            }
            if (idEntry.nonTight == "auauMLPSideband")
            {
                std::cout << " nonTight=(" << cfg.auau_nontight_mlp_min_slope << " * ET + " << cfg.auau_nontight_mlp_min_intercept
                          << ", " << cfg.auau_nontight_mlp_max_slope << " * ET + " << cfg.auau_nontight_mlp_max_intercept << ")";
            }
            else
            {
                std::cout << " nonTight=complement";
            }
            std::cout << "\n";
        }
        
        std::cout << "[CFG] isolation mode:";
        if (!entryIsSlidingIso)
        {
            std::cout << " fixed (thrReco=fixedGeV=" << entryFixedGeV << ")";
        }
        else
        {
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
            const auto& coneWps = (std::fabs(entryConeR - 0.30) < 0.015 && !cfg.auauCentIsoWPR30.empty())
                ? cfg.auauCentIsoWPR30
                : ((std::fabs(entryConeR - 0.40) < 0.015 && !cfg.auauCentIsoWPR40.empty())
                   ? cfg.auauCentIsoWPR40
                   : cfg.auauCentIsoWP);
            if (coneWps.size() == 1)
            {
                std::cout << " sliding -> AuAu/embedded centrality-fit"
                << " (thrReco(cent)=" << coneWps.front().aGeV
                << " + " << coneWps.front().bPerGeV << " * cent)";
            }
            else if (!coneWps.empty())
            {
                std::cout << " sliding -> per-centrality-bin WP list"
                << " (nCentWP=" << coneWps.size() << ", coneR=" << entryConeR << ")";
            }
            else
#endif
            {
                const auto* ppConeWP = (std::fabs(entryConeR - 0.30) < 0.015 && cfg.ppIsoWPR30.configured)
                    ? &cfg.ppIsoWPR30
                    : ((std::fabs(entryConeR - 0.40) < 0.015 && cfg.ppIsoWPR40.configured)
                       ? &cfg.ppIsoWPR40
                       : nullptr);
                if (ppConeWP)
                {
                    std::cout << " sliding -> pp cone-fit"
                    << " (coneR=" << entryConeR
                    << ", thrReco(pT)=" << ppConeWP->aGeV
                    << " + " << ppConeWP->bPerGeV << " * pTgamma"
                    << ", sideGap=" << ppConeWP->sideGapGeV << ")";
                }
                else
                {
                    std::cout << " sliding -> legacy global fallback"
                    << " (thrReco(pT)=" << cfg.isoA
                    << " + " << cfg.isoB << " * pTgamma)";
                }
            }
        }
        std::cout << " sideGap=" << cfg.isoGap << "\n";
    }
    
    recoilJets->enableEventDisplayDiagnostics(cfg.event_display_tree);
    recoilJets->setEventDisplayDiagnosticsMaxPerBin(cfg.event_display_tree_max_per_bin);
    
    if (vlevel > 0)
    {
        std::cout << "[CFG] EventDisplayTree: enable=" << (cfg.event_display_tree ? "true" : "false")
        << " max_per_bin=" << cfg.event_display_tree_max_per_bin << "\n";
    }
    
    recoilJets->Verbosity(vlevel);
    if (verbose) std::cout << "[INFO] RJ_VERBOSITY → " << vlevel << "\n";
    
    std::string dtype;
    if (isSim)
    {
        if (datasetToken == "issimembedded" || datasetToken == "simembedded" ||
            datasetToken == "issimembeddedinclusive" || datasetToken == "simembeddedinclusive")
            dtype = "isSimEmbedded";
        else if (datasetToken == "issimjet5" || datasetToken == "simjet5"
                 || datasetToken == "issiminclusive" || datasetToken == "siminclusive")
            dtype = "isSimInclusive";
        else if (datasetToken == "issimmb" || datasetToken == "simmb")
            dtype = "isSimMB";
        else
            dtype = "isSim";
    }
    else if (isAuAuData) dtype = "isAuAu";
    else dtype = "isPP";
    if (isPPrun25) dtype = "isPP";
#if defined(RJ_UNIFIED_ANALYSIS_PP)
    if (dtype == "isAuAu" || dtype == "isSimEmbedded")
        detail::bail("Fun4All_recoilJets.C wrapper is pp analysis-module only; use Fun4All_recoilJets_AuAu.C for isAuAu/isSimEmbedded jobs.");
#endif
    
    if (verbose || vlevel > 0)
    {
        std::cout << "[FLOW] final RecoilJets mode:"
        << " dtype=" << dtype
        << " | datasetToken=" << datasetToken
        << " | photonInputClusterNode=" << photonInputClusterNode
        << " | photonBuilderIsAuAu=" << (photonBuilderIsAuAu ? "true" : "false")
        << " | preselectionVariant=" << idEntry.preselection
        << " | tightVariant=" << idEntry.tight
        << " | preselectionPhotonNode=" << preselectionPhotonNode
        << " | tightPhotonNode=" << tightPhotonNode
        << std::endl;
    }
    recoilJets->setDataType(dtype);
    };

    if (usePPG12PhotonYieldTopoIso && !ppg12TopoBuilderRegistered)
    {
        registerPPG12PhotonYieldTopoBuilder("before RecoilJets fanout");
    }

    if (env_truthy_local("RJ_PPG12_FIG8_BUILD_NOSPLIT") &&
        env_truthy_local("RJ_PPG12_PHOTON_YIELD") &&
        isSim && !isAuAuLike &&
        !ppg12ClosureCanary)
    {
        auto* ppg12Fig8NoSplitBuilder =
            new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate_PPG12Fig8NoSplit");
        ppg12Fig8NoSplitBuilder->Detector("CEMC");
        ppg12Fig8NoSplitBuilder->set_threshold_energy(0.070);
        const char* calibroot = std::getenv("CALIBRATIONROOT");
        if (calibroot && std::string(calibroot).size())
        {
            const std::string emcProfile =
                std::string(calibroot) + "/EmcProfile/CEMCprof_Thresh30MeV.root";
            ppg12Fig8NoSplitBuilder->LoadProfile(emcProfile);
        }
        ppg12Fig8NoSplitBuilder->setSubclusterSplitting(false);
        ppg12Fig8NoSplitBuilder->setOutputClusterNodeName("CLUSTERINFO_CEMC_NO_SPLIT");
        ppg12Fig8NoSplitBuilder->set_UseTowerInfo(1);
        ppg12Fig8NoSplitBuilder->Verbosity(0);
        se->registerSubsystem(ppg12Fig8NoSplitBuilder);
        if (vlevel > 0)
        {
            std::cout << "[PPG12_FIG8] registered no-split template cluster builder"
                      << " node=CLUSTERINFO_CEMC_NO_SPLIT"
                      << " threshold=0.070"
                      << " profile=CEMCprof_Thresh30MeV.root"
                      << " before RecoilJets fanout\n";
        }
    }
    
    std::vector<RecoilJets*> recoilJetsInstances;
    recoilJetsInstances.reserve(idFanoutEntries.size());
    for (std::size_t i = 0; i < idFanoutEntries.size(); ++i)
    {
        const auto& entry = idFanoutEntries[i];
        std::ostringstream moduleName;
        moduleName << "RecoilJets_ID" << i;
        auto* recoilJets = new RecoilJets(entry.outRoot, moduleName.str());
        configureRecoilJetsInstance(recoilJets, entry);
        se->registerSubsystem(recoilJets);
        recoilJetsInstances.push_back(recoilJets);
        if (vlevel > 0)
        {
            std::cout << "[ID-FANOUT] registered " << moduleName.str()
                      << " cfgTag=" << entry.cfgTag
                      << " output=" << entry.outRoot
                      << " preselection=" << entry.preselection
                      << " tight=" << entry.tight
                      << " nonTight=" << entry.nonTight << "\n";
        }
    }
    
    //--------------------------------------------------------------------
    // 6.  Run
    //--------------------------------------------------------------------
    try
    {
        if (vlevel > 0) std::cout << "[INFO] Starting event loop …" << std::endl;
        
        const bool stepEvents = ([]{
            const char* env = std::getenv("RJ_STEP_EVENTS");
            return (env && std::atoi(env) != 0);
        })();
        int runRc = Fun4AllReturnCodes::EVENT_OK;
        
        if (stepEvents && nEvents > 0)
        {
            for (int ievt = 0; ievt < nEvents; ++ievt)
            {
                std::cout << "[RUN] >>> event " << (ievt + 1) << "/" << nEvents << std::endl;
                if (const char* trace = std::getenv("RJ_REPLAY_TRACE"); trace && std::atoi(trace) != 0)
                {
                    se->Verbosity(20);
                    std::cout << "[RUN] subsystem tracing enabled for event " << (ievt + 1) << std::endl;
                }
                const int rc = se->run(1);
                runRc = rc;
                std::cout << "[RUN] <<< event " << (ievt + 1) << "/" << nEvents << "  rc=" << rc << std::endl;
                if (rc != 0) break;
            }
        }
        else
        {
            runRc = se->run(nEvents);
        }
        
        // RecoilJets AuAu centrality counters do not exist in the pp class.
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
        std::uint64_t centralityValidTotal = 0;
        std::uint64_t centralityInvalidSkippedTotal = 0;
        std::uint64_t centralityStaleGuardedTotal = 0;
        if (isAuAuData)
        {
            for (std::size_t i = 0; i < recoilJetsInstances.size(); ++i)
            {
                const auto* recoilJets = recoilJetsInstances[i];
                const auto valid = recoilJets->validCentralityObservedEvents();
                const auto invalid = recoilJets->invalidCentralityObservedEvents();
                const auto staleGuarded = recoilJets->staleCentralityGuardedEvents();
                const auto evaluated = valid + invalid;
                centralityValidTotal += valid;
                centralityInvalidSkippedTotal += invalid;
                centralityStaleGuardedTotal += staleGuarded;
                const char* status =
                    (valid > 0 && invalid > 0) ? "PASS_WITH_SKIPPED_INVALID" :
                    (valid > 0) ? "PASS" :
                    (invalid > 0) ? "EMPTY_ALL_INVALID_EVENTS_SKIPPED" :
                    "EMPTY_NO_EVALUABLE_EVENTS";
                std::cout << "[AUAU_CENTRALITY_CONTRACT] module=RecoilJets_ID" << i
                          << " valid=" << valid
                          << " invalid=" << invalid
                          << " evaluated=" << evaluated
                          << " stale_guarded=" << staleGuarded
                          << " status=" << status << std::endl;
            }
            std::cout << "[AUAU_CENTRALITY_CONTRACT_SUMMARY] valid=" << centralityValidTotal
                      << " invalid_skipped=" << centralityInvalidSkippedTotal
                      << " stale_guarded=" << centralityStaleGuardedTotal
                      << " action=invalid_events_audited_and_skipped" << std::endl;
        }
#endif

        // Keep the event loop quiet in batch, but surface the bounded terminal
        // subsystem diagnostics before interpreting Fun4All's aggregate End
        // return. This does not change event processing or ROOT content.
        _silence.disable();
        if (vlevel > 0) std::cout << "[INFO] Calling se->End() …" << std::endl;
        const int endRc = se->End();
        detail::enforce_fun4all_status(
            "analysis",
            se,
            nEvents,
            runRc,
            endRc,
            permittedRepeatingPedestalInputManager);
        if (vlevel > 0) std::cout << "[INFO] Finished successfully." << std::endl;
    }
    catch (const std::exception& e)
    {
        detail::bail(std::string("exception in Fun4All: ") + e.what());
    }
}

#endif   // ROOT_VERSION guard
