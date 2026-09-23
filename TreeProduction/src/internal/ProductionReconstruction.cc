// ProductionReconstruction.cc
//
// Which modules build the event before PhotonJetTree sees it?
// Owns the complete upstream reconstruction graph in physical registration order.
// Does not own configuration parsing or physical input streams.
// Called through Production.h before PhotonJetTree captures the event.
//
#include "Production.h"
#include "../PhotonJetTree.h"

// ---- sPHENIX framework -----------------------------------------------------

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/CdbUrlSave.h>
#include <ffaobjects/EventHeader.h>
#include <ffaobjects/RunHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <cdbobjects/CDBTTree.h>

// ---- calorimeter -------------------------------------------------------------

#include <calobase/RawTowerDefs.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/PhotonClusterBuilder.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterBuilderTopo.h>

#include <calostatusskimmer/CaloStatusSkimmer.h>

// ---- event geometry, minimum bias, centrality -------------------------------

#include <calotrigger/MinimumBiasClassifier.h>
#include <centrality/CentralityReco.h>
#include <g4mbd/MbdDigitization.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexReco.h>
#include <mbd/MbdReco.h>
#include <zdcinfo/ZdcReco.h>

// ---- jets ---------------------------------------------------------------------

#include <g4jets/TruthJetInput.h>

#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>

#include <jetbase/FastJetOptions.h>
#include <jetbase/Jet.h>
#include <jetbase/JetAlgo.h>
#include <jetbase/JetCalib.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

// ---- ROOT and standard library ------------------------------------------------

#include <TFile.h>
#include <phool/RunnumberRange.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
using photonjet::production::Conditions;
using photonjet::production::InputLayout;
using photonjet::production::InputStream;
using photonjet::production::Job;
using photonjet::production::Plan;
using photonjet::production::Reconstruction;
using photonjet::production::TruthJetMode;
using photonjet::production::detail::fail;
using photonjet::production::detail::requireStatus;
using photonjet::production::detail::kCemcChannelCount;
using photonjet::production::detail::kArchivedDoubleInteractionRun;
using photonjet::production::detail::kDataRunThreshold;

constexpr int kArchivedDoubleInteractionPedestalSequence = 534;

constexpr std::array<unsigned int, 5> kArchivedDoubleInteractionSeeds =
{
  2991264730U, 4256268992U, 2394322166U, 874466025U, 2240380304U
};

// Exposure observer
//
// Registered before every other module. It reports each input event to the
// producer before reconstruction has a chance to abort it, which is the only
// way the producer can count the events it never sees.
//
class ExposureObserver final : public SubsysReco
{
 public:
  ExposureObserver(PhotonJetTree* producer, std::string eventHeaderNode)
    : SubsysReco("PhotonJetTreeExposureObserver")
    , m_producer(producer)
    , m_eventHeaderNode(std::move(eventHeaderNode))
  {
  }

  int process_event(PHCompositeNode* topNode) override
  {
    int run = 0;
    std::int64_t sequence = -1;

    if (auto* header = findNode::getClass<EventHeader>(topNode, m_eventHeaderNode);
        header && header->isValid())
    {
      run = header->get_RunNumber();
      sequence = header->get_EvtSequence();
    }

    m_producer->observeUpstreamEvent(run, sequence);
    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  PhotonJetTree* m_producer;
  std::string m_eventHeaderNode;
};

//
// Jet clustering with active area
//
// The release FastJetAlgoSub keeps negative-energy tower constituents by
// clustering them as tiny positive pseudojets and restoring the original
// four-vectors afterwards, but it does not compute the FastJet active area.
// This class keeps that clustering contract exactly and adds the active-area
// ghosts and the standard area property, which the stored jet area requires.
//
// TODO(upstream): propose area support for FastJetAlgoSub in coresoftware so
// this class can be retired; until then it stays with the producer that
// needs it.
//
class JetAreaAlgorithm final : public JetAlgo
{
 public:
  explicit JetAreaAlgorithm(const FastJetOptions& options)
    : m_options(options)
  {
    // Print the FastJet banner once, silently unless verbose.
    if (m_options.verbosity > 0)
    {
      fastjet::ClusterSequence::print_banner();
    }
    else
    {
      std::ostringstream sink;
      fastjet::ClusterSequence::set_fastjet_banner_stream(&sink);
      fastjet::ClusterSequence::print_banner();
      fastjet::ClusterSequence::set_fastjet_banner_stream(&std::cout);
    }
  }

  ~JetAreaAlgorithm() override = default;

  void identify(std::ostream& os = std::cout) override
  {
    os << "   JetAreaAlgorithm: ";
    if (m_options.algo == Jet::ANTIKT) os << "ANTIKT";
    else if (m_options.algo == Jet::KT) os << "KT";
    else if (m_options.algo == Jet::CAMBRIDGE) os << "CAMBRIDGE";
    else os << "UNKNOWN";
    os << " r=" << m_options.jet_R
       << " active_area=" << (m_options.calc_area ? 1 : 0) << std::endl;
  }

  Jet::ALGO get_algo() override { return m_options.algo; }
  float get_par() override { return m_options.jet_R; }

  void cluster_and_fill(std::vector<Jet*>& particles, JetContainer* jetcont) override
  {
    if (!jetcont)
    {
      return;
    }
    initializeContainer(jetcont);

    std::vector<fastjet::PseudoJet> pseudojets;
    pseudojets.reserve(particles.size());
    for (std::size_t ipart = 0; ipart < particles.size(); ++ipart)
    {
      Jet* particle = particles[ipart];
      if (!particle)
      {
        continue;
      }
      float energy = particle->get_e();
      if (energy == 0.0F)
      {
        continue;
      }
      float px = particle->get_px();
      float py = particle->get_py();
      float pz = particle->get_pz();

      // The negative-energy contract of the release algorithm.
      if (energy < 0.0F)
      {
        const float ratio = 0.001F / energy;
        energy *= ratio;
        px *= ratio;
        py *= ratio;
        pz *= ratio;
      }

      fastjet::PseudoJet pseudojet(px, py, pz, energy);
      pseudojet.set_user_index(static_cast<int>(ipart));
      pseudojets.push_back(pseudojet);
    }

    const fastjet::JetDefinition definition = jetDefinition();
    std::unique_ptr<fastjet::ClusterSequence> plainSequence;
    std::unique_ptr<fastjet::ClusterSequenceArea> areaSequence;
    std::vector<fastjet::PseudoJet> fastjets;

    if (m_options.calc_area)
    {
      const fastjet::AreaDefinition areaDefinition(
          fastjet::active_area_explicit_ghosts,
          fastjet::GhostedAreaSpec(m_options.ghost_max_rap, 1, m_options.ghost_area));
      areaSequence = std::make_unique<fastjet::ClusterSequenceArea>(
          pseudojets, definition, areaDefinition);
      fastjets = areaSequence->inclusive_jets();
    }
    else
    {
      plainSequence = std::make_unique<fastjet::ClusterSequence>(pseudojets, definition);
      fastjets = plainSequence->inclusive_jets();
    }

    unsigned int outputId = 0;
    for (const fastjet::PseudoJet& fastjet : fastjets)
    {
      if (m_options.calc_area && fastjet.is_pure_ghost())
      {
        continue;
      }

      Jet* jet = jetcont->add_jet();
      if (!jet)
      {
        continue;
      }

      float totalPx = 0.0F, totalPy = 0.0F, totalPz = 0.0F, totalEnergy = 0.0F;
      float weightedTime = 0.0F, weightedEnergy = 0.0F;

      for (const fastjet::PseudoJet& constituent : fastjet.constituents())
      {
        if (m_options.calc_area && constituent.is_pure_ghost())
        {
          continue;
        }
        const int inputIndex = constituent.user_index();
        if (inputIndex < 0 || inputIndex >= static_cast<int>(particles.size()))
        {
          continue;
        }
        Jet* particle = particles[static_cast<std::size_t>(inputIndex)];
        if (!particle)
        {
          continue;
        }

        totalPx += particle->get_px();
        totalPy += particle->get_py();
        totalPz += particle->get_pz();
        totalEnergy += particle->get_e();
        if (particle->size_properties() > Jet::PROPERTY::prop_t &&
            !std::isnan(particle->get_property(Jet::PROPERTY::prop_t)))
        {
          weightedTime += particle->get_property(Jet::PROPERTY::prop_t) * particle->get_e();
          weightedEnergy += particle->get_e();
        }
        jet->insert_comp(particle->get_comp_vec(), true);
      }

      if (jet->size_properties() < Jet::PROPERTY::prop_t + 1)
      {
        jet->resize_properties(Jet::PROPERTY::prop_t + 1);
      }
      jet->set_property(Jet::PROPERTY::prop_t, weightedTime / weightedEnergy);
      if (m_options.calc_area)
      {
        jet->set_property(m_areaIndex, fastjet.area());
      }

      jet->set_comp_sort_flag();
      jet->set_px(totalPx);
      jet->set_py(totalPy);
      jet->set_pz(totalPz);
      jet->set_e(totalEnergy);
      jet->set_id(outputId++);
    }
  }

 private:
  void initializeContainer(JetContainer* jetcont)
  {
    if (m_initialized)
    {
      return;
    }
    m_options.initialize();
    jetcont->set_algo(m_options.algo);
    jetcont->set_jetpar_R(m_options.jet_R);
    if (m_options.calc_area)
    {
      jetcont->add_property(Jet::PROPERTY::prop_area);
      m_areaIndex = jetcont->property_index(Jet::PROPERTY::prop_area);
    }
    m_initialized = true;
  }

  fastjet::JetDefinition jetDefinition() const
  {
    if (m_options.algo == Jet::KT)
    {
      return fastjet::JetDefinition(fastjet::kt_algorithm, m_options.jet_R,
                                    fastjet::E_scheme, fastjet::Best);
    }
    if (m_options.algo == Jet::CAMBRIDGE)
    {
      return fastjet::JetDefinition(fastjet::cambridge_algorithm, m_options.jet_R,
                                    fastjet::E_scheme, fastjet::Best);
    }
    return fastjet::JetDefinition(fastjet::antikt_algorithm, m_options.jet_R,
                                  fastjet::E_scheme, fastjet::Best);
  }

  FastJetOptions m_options{};
  bool m_initialized = false;
  Jet::PROPERTY m_areaIndex = Jet::PROPERTY::no_property;
};

JetAreaAlgorithm* antiKt(const float radius, const bool withArea)
{
  FastJetOptions options{};
  options.algo = Jet::ANTIKT;
  options.jet_R = radius;
  options.use_jet_min_pt = true;
  options.jet_min_pt = 0.0F;
  options.calc_area = withArea;
  options.verbosity = 0;
  return new JetAreaAlgorithm(options);
}

//
// Jet calibration node guard
//
// JetCalib expects the DST/ANTIKT and DST/TOWER composite nodes to exist.
// Some inputs do not carry them. This creates only the missing containers and
// touches no physics object.
//
class EnsureJetCalibNodes final : public SubsysReco
{
 public:
  EnsureJetCalibNodes() : SubsysReco("PhotonJetTreeEnsureJetCalibNodes") {}

  int InitRun(PHCompositeNode* topNode) override
  {
    PHNodeIterator iterator(topNode);
    auto* dst = dynamic_cast<PHCompositeNode*>(iterator.findFirst("PHCompositeNode", "DST"));
    if (!dst)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }
    for (const char* name : {"ANTIKT", "TOWER"})
    {
      PHNodeIterator local(dst);
      if (!dynamic_cast<PHCompositeNode*>(local.findFirst("PHCompositeNode", name)))
      {
        dst->addNode(new PHCompositeNode(name));
      }
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }
};

//
// CEMC tower-status contract, Au+Au data
//
// Paired Au+Au data files carry the raw tower status that the calibration
// copies into the calibrated towers. Two guards bracket the calibration: the
// first validates the raw status against the bad-tower map the original file
// was reconstructed with, restoring missing hot flags only when the packet
// explicitly authorises a run-bound recovery payload; the second confirms the
// calibrated towers carry the same rejection. Neither guard reads the current
// conditions tag: the map comes from the original file's saved provenance or
// from the pinned payload, never from a lookup that could differ from what
// the reconstruction actually used.
//
constexpr std::uint8_t kCemcHotBit = 1;

struct CemcStatusState
{
  int expectedRun = 0;
  std::string originalCalorimeterInput;
  bool recoveryAuthorized = false;
  std::string recoveryPayload;
  std::string selectedPayload;
  bool restoringMissingStatus = false;
  bool initialized = false;
  std::vector<std::uint8_t> mask;
  unsigned long long rawEvents = 0;
  unsigned long long calibratedEvents = 0;
};

class CemcTowerStatusGuard final : public SubsysReco
{
 public:
  CemcTowerStatusGuard(std::shared_ptr<CemcStatusState> state, const bool calibrated)
    : SubsysReco(calibrated ? "CemcStatusBeforeClustering" : "CemcStatusBeforeCalibration")
    , m_state(std::move(state))
    , m_calibrated(calibrated)
  {
  }

  int InitRun(PHCompositeNode* top) override
  {
    try
    {
      requireStatus(bool(m_state) && m_state->expectedRun > 0, "missing run binding");
      auto* runHeader = findNode::getClass<RunHeader>(top, "RunHeader");
      requireStatus(runHeader && runHeader->get_RunNumber() == m_state->expectedRun,
                    "input RunHeader does not match the bound run");
      if (m_calibrated)
      {
        requireStatus(m_state->initialized, "calibrated check registered before raw check");
        return Fun4AllReturnCodes::EVENT_OK;
      }

      auto* towers = findNode::getClass<TowerInfoContainer>(top, "TOWERS_CEMC");
      requireStatus(towers && towers->size() == kCemcChannelCount, "TOWERS_CEMC input missing");

      // Read only the original calorimeter file's RUN metadata: the live RUN
      // node is shared with the paired and geometry inputs and its last
      // loaded conditions URL is not necessarily the calorimeter's.
      requireStatus(!m_state->originalCalorimeterInput.empty(),
                    "missing exact original calorimeter source");
      auto originalRun = std::make_unique<PHCompositeNode>("ORIGINAL_CALO_RUN");
      {
        PHNodeIOManager original(m_state->originalCalorimeterInput, PHReadOnly, PHRunTree);
        requireStatus(original.isFunctional() && original.read(originalRun.get()),
                      "unreadable original calorimeter RUN metadata");
      }
      auto* originalHeader = findNode::getClass<RunHeader>(originalRun.get(), "RunHeader");
      requireStatus(originalHeader && originalHeader->get_RunNumber() == m_state->expectedRun,
                    "original calorimeter metadata run mismatch");

      std::set<std::string> saved;
      if (auto* urls = findNode::getClass<CdbUrlSave>(originalRun.get(), "CdbUrl"))
      {
        for (auto it = urls->begin(); it != urls->end(); ++it)
        {
          if (std::get<0>(*it) == "CEMC_BadTowerMap")
          {
            saved.insert(std::get<1>(*it));
          }
        }
      }
      requireStatus(saved.size() <= 1, "ambiguous saved CEMC map provenance");

      if (saved.empty())
      {
        requireStatus(m_state->recoveryAuthorized && !m_state->recoveryPayload.empty(),
                      "no saved CEMC map; explicit validated recovery payload required");
        // The repair covers the demonstrated missing-map, zero-hot-flag case only.
        for (unsigned i = 0; i < towers->size(); ++i)
        {
          auto* t = towers->get_tower_at_channel(i);
          requireStatus(t && !t->get_isHot(),
                        "missing provenance with existing HOT state requires review");
        }
        m_state->selectedPayload = m_state->recoveryPayload;
        m_state->restoringMissingStatus = true;
      }
      else
      {
        // An explicit repair payload never refreshes already valid inputs.
        m_state->selectedPayload = *saved.begin();
        m_state->restoringMissingStatus = false;
      }

      const std::string& path = m_state->selectedPayload;
      const std::string suffix = "_" + std::to_string(m_state->expectedRun) + "cdb.root";
      requireStatus(path.find("/cdb/CEMC_BadTowerMap/") != std::string::npos &&
                        path.size() > suffix.size() &&
                        path.compare(path.size() - suffix.size(), suffix.size(), suffix) == 0,
                    "payload path is not bound to this run and the CEMC domain");

      std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
      requireStatus(file && !file->IsZombie(), "unreadable CEMC map file");
      file->Close();

      CDBTTree map(path);
      m_state->mask.clear();
      m_state->mask.reserve(kCemcChannelCount);
      for (unsigned i = 0; i < towers->size(); ++i)
      {
        const int value = map.GetIntValue(towers->encode_key(i), "status", 0);
        // A missing integer field reads as INT_MIN and must never count as good.
        requireStatus(value >= 0 && value <= 255, "missing or invalid map status field");
        m_state->mask.push_back(value > 0);
      }
      m_state->initialized = true;

      std::cout << "CEMC_STATUS_CONTRACT_INIT run=" << m_state->expectedRun
                << " original_input=" << m_state->originalCalorimeterInput
                << " payload=" << path
                << " action=" << (m_state->restoringMissingStatus ? "RESTORE_MISSING_HOT_ONLY"
                                                                  : "PRESERVE_VALIDATE")
                << " map_rejected=" << std::count(m_state->mask.begin(), m_state->mask.end(), 1)
                << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
    catch (const std::exception& error)
    {
      return fatal(error.what());
    }
  }

  int process_event(PHCompositeNode* top) override
  {
    try
    {
      requireStatus(m_state && m_state->initialized, "uninitialized status guard");
      auto* towers = findNode::getClass<TowerInfoContainer>(
          top, m_calibrated ? "TOWERINFO_CALIB_CEMC" : "TOWERS_CEMC");
      requireStatus(towers && towers->size() == kCemcChannelCount, "required CEMC node missing");

      if (!m_calibrated && m_state->restoringMissingStatus)
      {
        for (unsigned i = 0; i < towers->size(); ++i)
        {
          if (m_state->mask[i])
          {
            towers->get_tower_at_channel(i)->set_isHot(true);
          }
        }
      }

      for (unsigned i = 0; i < towers->size(); ++i)
      {
        auto* t = towers->get_tower_at_channel(i);
        requireStatus(t != nullptr, "null tower");
        requireStatus(!m_state->mask[i] || (t->get_status() & kCemcHotBit),
                      "map-rejected tower lacks HOT flag at channel " + std::to_string(i));
        requireStatus(!m_state->mask[i] || (t->get_isHot() && !t->get_isGood()),
                      "map-rejected tower must be HOT and not GOOD in the actual node");
      }

      if (m_calibrated)
      {
        ++m_state->calibratedEvents;
        requireStatus(m_state->calibratedEvents == m_state->rawEvents,
                      "raw and calibrated status event ordering mismatch");
      }
      else
      {
        ++m_state->rawEvents;
      }
      return Fun4AllReturnCodes::EVENT_OK;
    }
    catch (const std::exception& error)
    {
      return fatal(error.what());
    }
  }

  int End(PHCompositeNode*) override
  {
    if (m_calibrated)
    {
      if (!m_state || !m_state->initialized || m_state->rawEvents != m_state->calibratedEvents)
      {
        return fatal("incomplete paired CEMC status audit");
      }
      std::cout << "CEMC_STATUS_CONTRACT_END run=" << m_state->expectedRun
                << " raw_events=" << m_state->rawEvents
                << " calibrated_events=" << m_state->calibratedEvents << " PASS" << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  static int fatal(const std::string& reason)
  {
    std::cerr << "CEMC_STATUS_CONTRACT_FATAL " << reason << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::shared_ptr<CemcStatusState> m_state;
  bool m_calibrated;
};
//
// Tower jet inputs
//
void addTowerInputs(JetReco* reco, const photonjet::JetView view, const bool auau)
{
  Jet::SRC cemc = Jet::CEMC_TOWERINFO;
  Jet::SRC ihcal = Jet::HCALIN_TOWERINFO;
  Jet::SRC ohcal = Jet::HCALOUT_TOWERINFO;

  if (view == photonjet::JetView::AuAuSub1)
  {
    cemc = Jet::CEMC_TOWERINFO_SUB1;
    ihcal = Jet::HCALIN_TOWERINFO_SUB1;
    ohcal = Jet::HCALOUT_TOWERINFO_SUB1;
  }
  else if (view == photonjet::JetView::AuAuNoSub)
  {
    cemc = Jet::CEMC_TOWERINFO_RETOWER;
  }

  for (const Jet::SRC source : {cemc, ihcal, ohcal})
  {
    auto* input = new TowerJetInput(source, "TOWERINFO_CALIB");
    if (auau)
    {
      input->set_GlobalVertexType(GlobalVertex::MBD);
    }
    reco->add_input(input);
  }
}

//
// Registration, one function per reconstruction stage
//
// Minimum-bias detector, zero-degree calorimeter and global vertex.
void registerEventGeometry(Fun4AllServer* server, const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  if (!reco.reconstructEventGeometry)
  {
    return;
  }

  if (reco.isArchivedDoubleInteraction)
  {
    // The archived Geant4-only input has no digitised minimum-bias detector.
    server->registerSubsystem(new MbdDigitization());
  }

  server->registerSubsystem(new MbdReco());


  server->registerSubsystem(new GlobalVertexReco());
}

// The archived double-interaction reconstruction replays a fixed random
// sequence so that the reconstructed calorimeter matches the donor exactly.
void registerArchivedDoubleInteractionSeeds(const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  if (!reco.isArchivedDoubleInteraction)
  {
    return;
  }

  recoConsts* rc = recoConsts::instance();
  if (reco.forbidRecoConstsRandomSeed && rc->FlagExist("RANDOMSEED"))
  {
    fail("archived double-interaction reconstruction forbids recoConsts RANDOMSEED");
  }
  if (reco.historicalSeeds.size() != kArchivedDoubleInteractionSeeds.size())
  {
    fail("archived double-interaction reconstruction requires the exact five-seed sequence");
  }
  for (std::size_t i = 0; i < kArchivedDoubleInteractionSeeds.size(); ++i)
  {
    if (reco.historicalSeeds[i] != kArchivedDoubleInteractionSeeds[i])
    {
      fail("archived double-interaction seed sequence differs from the frozen sequence");
    }
  }
  if (reco.expectedPedestalSequence != kArchivedDoubleInteractionPedestalSequence)
  {
    fail("archived double-interaction pedestal sequence must be 534");
  }

  PHRandomSeed::Verbosity(1);
  for (const unsigned int seed : reco.historicalSeeds)
  {
    PHRandomSeed::LoadSeed(seed);
  }
}

// Calorimeter calibration, status and clustering, in the order of the
// historical Process_Calo_Calib.
void registerCalorimeter(Fun4AllServer* server, const Plan& plan,
                         const std::shared_ptr<CemcStatusState>& cemcStatus)
{
  const Reconstruction& reco = plan.reconstruction;
  const Conditions& cond = plan.conditions;
  if (!reco.reconstructCalorimeter)
  {
    return;
  }

  const int timestamp = static_cast<int>(cond.cdbTimestamp);
  const bool simulationTimestamp = timestamp <= kDataRunThreshold;

  // Drop events the event combiner left incomplete. Data only.
  if (reco.runCalorimeterStatusSkimmer)
  {
    server->registerSubsystem(new CaloStatusSkimmer("CaloStatusSkimmer"));
  }

  // Zero-degree calorimeter towers, needed by the minimum-bias classifier.
  if (reco.reconstructZdc)
  {
    auto* zdc = new CaloTowerBuilder("ZDCBUILDER");
    zdc->set_detector_type(CaloTowerDefs::ZDC);
    zdc->set_builder_type(CaloTowerDefs::kPRDFTowerv4);
    const bool fastProcessing =
        (timestamp > RunnumberRange::RUN2PP_FIRST && timestamp < RunnumberRange::RUN2PP_LAST) ||
        (timestamp > RunnumberRange::RUN3PP_FIRST && timestamp < RunnumberRange::RUN3PP_LAST);
    if (fastProcessing)
    {
      zdc->set_processing_type(CaloWaveformProcessing::FAST);
    }
    else
    {
      zdc->set_processing_type(CaloWaveformProcessing::FUNCFIT);
      zdc->set_funcfit_type(2);
    }
    zdc->set_nsamples(16);
    zdc->set_offlineflag();
    server->registerSubsystem(zdc);
  }

  if (reco.reconstructZdc)
  {
    auto* zdc = new ZdcReco();
    zdc->set_zdc1_cut(0.0);
    zdc->set_zdc2_cut(0.0);
    server->registerSubsystem(zdc);
  }

  // Raw-tower status validation against the original bad-tower map.
  if (cemcStatus)
  {
    server->registerSubsystem(new CemcTowerStatusGuard(cemcStatus, false));
  }

  // Legacy status setters. Embedded simulation keeps its producer status.
  if (reco.setTowerStatus)
  {
    auto* statusEmc = new CaloTowerStatus("CEMCSTATUS");
    statusEmc->set_detector_type(CaloTowerDefs::CEMC);
    if (simulationTimestamp)
    {
      const std::string hotMap = CDBInterface::instance()->getUrl(reco.simulationHotTowerMapKey);
      statusEmc->set_directURL_hotMap(hotMap);
    }
    server->registerSubsystem(statusEmc);

    auto* statusIhcal = new CaloTowerStatus("HCALINSTATUS");
    statusIhcal->set_detector_type(CaloTowerDefs::HCALIN);
    server->registerSubsystem(statusIhcal);

    auto* statusOhcal = new CaloTowerStatus("HCALOUTSTATUS");
    statusOhcal->set_detector_type(CaloTowerDefs::HCALOUT);
    server->registerSubsystem(statusOhcal);
  }

  // Tower calibration. The CEMC calibrator pins the zero-suppression cross
  // calibration payload when the profile asks for it.
  auto* calibEmc = new CaloTowerCalib("CEMCCALIB");
  calibEmc->set_detector_type(CaloTowerDefs::CEMC);
  if (cond.cemcZeroSuppressionCrossCalibration)
  {
    calibEmc->set_doZScrosscalib(true);
    calibEmc->set_doAbortNoZSCalib(true);
    calibEmc->set_directURL_ZScrosscalib(cond.cemcZeroSuppressionPayload);
  }
  server->registerSubsystem(calibEmc);

  auto* calibOhcal = new CaloTowerCalib("HCALOUT");
  calibOhcal->set_detector_type(CaloTowerDefs::HCALOUT);
  server->registerSubsystem(calibOhcal);

  auto* calibIhcal = new CaloTowerCalib("HCALIN");
  calibIhcal->set_detector_type(CaloTowerDefs::HCALIN);
  server->registerSubsystem(calibIhcal);

  if (cemcStatus)
  {
    server->registerSubsystem(new CemcTowerStatusGuard(cemcStatus, true));
  }

  // Simulation before run 28 carried a separate CEMC recalibration. Later
  // simulation moved it into the waveform stage.
  if (simulationTimestamp && timestamp < kArchivedDoubleInteractionRun)
  {
    const std::string recalibration = CDBInterface::instance()->getUrl("CEMC_MC_RECALIB");
    if (recalibration.empty())
    {
      fail("no CEMC simulation recalibration payload for timestamp " + std::to_string(timestamp));
    }
    auto* calibEmcMc = new CaloTowerCalib("CEMCCALIB_MC");
    calibEmcMc->set_detector_type(CaloTowerDefs::CEMC);
    calibEmcMc->set_inputNodePrefix("TOWERINFO_CALIB_");
    calibEmcMc->set_outputNodePrefix("TOWERINFO_CALIB_");
    calibEmcMc->set_directURL(recalibration);
    calibEmcMc->set_doCalibOnly(true);
    server->registerSubsystem(calibEmcMc);
  }

  // CEMC clusters.
  auto* clusters = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  clusters->Detector("CEMC");
  clusters->set_threshold_energy(reco.clusterThresholdGeV);
  clusters->LoadProfile(reco.clusterProfile);
  clusters->set_UseTowerInfo(1);
  clusters->set_UseAltZVertex(1);
  server->registerSubsystem(clusters);
}

// Au+Au data: minimum-bias classification and centrality from run-bound
// frozen payloads. Embedded simulation keeps the centrality of its input.
void registerCentrality(Fun4AllServer* server, const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  if (!reco.reconstructCentrality)
  {
    return;
  }

  if (reco.runMinimumBiasClassifier)
  {
    auto* minimumBias = new MinimumBiasClassifier();
    minimumBias->Verbosity(0);
    server->registerSubsystem(minimumBias);
  }

  auto* centrality = new CentralityReco();
  centrality->Verbosity(0);
  if (plan.conditions.centralityLocalPayloads)
  {
    centrality->setOverwriteDivs(plan.conditions.centralityDivisions);
    centrality->setOverwriteScale(plan.conditions.centralityRunScale);
    centrality->setOverwriteVtx(plan.conditions.centralityVertexScale);
  }
  server->registerSubsystem(centrality);
}

// Au+Au underlying-event subtraction: retower the CEMC onto the hadronic
// segmentation, estimate the background from seed jets in two passes, and
// write the subtracted tower containers the Sub1 jets are built from.
void registerAuAuBackground(Fun4AllServer* server, const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  if (!reco.reconstructAuAuBackground)
  {
    return;
  }

  auto* retower = new RetowerCEMC();
  retower->set_towerinfo(true);
  retower->set_frac_cut(reco.retowerFractionCut);
  retower->set_towerNodePrefix("TOWERINFO_CALIB");
  server->registerSubsystem(retower);

  auto* seedJets = new JetReco("JetsReco_HIRecoSeedsRaw_r02");
  addTowerInputs(seedJets, photonjet::JetView::AuAuNoSub, true);
  seedJets->add_algo(antiKt(0.2F, false), "AntiKt_TowerInfo_HIRecoSeedsRaw_r02");
  seedJets->set_algo_node("ANTIKT");
  seedJets->set_input_node("TOWER");
  server->registerSubsystem(seedJets);

  auto* background1 = new DetermineTowerBackground();
  background1->SetBackgroundOutputName("TowerInfoBackground_Sub1");
  background1->SetFlow(reco.flowModulation);
  background1->SetSeedType(0);
  background1->SetSeedJetD(reco.firstPassSeedJetD);
  background1->set_towerNodePrefix("TOWERINFO_CALIB");
  server->registerSubsystem(background1);

  auto* copySubtract = new CopyAndSubtractJets();
  copySubtract->SetFlowModulation(reco.flowModulation);
  copySubtract->set_towerinfo(true);
  copySubtract->set_towerNodePrefix("TOWERINFO_CALIB");
  server->registerSubsystem(copySubtract);

  auto* background2 = new DetermineTowerBackground();
  background2->SetBackgroundOutputName("TowerInfoBackground_Sub2");
  background2->SetFlow(reco.flowModulation);
  background2->SetSeedType(1);
  background2->SetSeedJetPt(reco.secondPassSeedPtGeV);
  background2->set_towerNodePrefix("TOWERINFO_CALIB");
  server->registerSubsystem(background2);

  auto* subtract = new SubtractTowers();
  subtract->SetFlowModulation(reco.flowModulation);
  subtract->set_towerinfo(true);
  subtract->set_towerNodePrefix("TOWERINFO_CALIB");
  server->registerSubsystem(subtract);
}

// Reconstructed jets, one raw collection per view and radius, each followed
// by exactly one jet energy-scale application into its corrected node.
void registerRecoJets(Fun4AllServer* server, const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  if (!reco.reconstructJets)
  {
    return;
  }

  server->registerSubsystem(new EnsureJetCalibNodes());

  for (const photonjet::JetNodeConfig& collection : plan.tree.jets.nodes)
  {
    std::ostringstream suffix;
    suffix << collection.inputIdentity << "_R" << std::lround(collection.radius * 100.0);

    auto* jetReco = new JetReco("JetReco_" + suffix.str());
    addTowerInputs(jetReco, collection.view, reco.isAuAu);
    jetReco->add_algo(antiKt(static_cast<float>(collection.radius), true), collection.rawNode);
    jetReco->set_algo_node("ANTIKT");
    jetReco->set_input_node(reco.isAuAu ? "TOWER" : "TOWERINFO_CALIB");
    server->registerSubsystem(jetReco);

    if (!plan.conditions.applyJetEnergyScale)
    {
      continue;
    }

    // The legacy method is selected explicitly. The payload it consumes was
    // resolved and hash-verified in ConfigureConditions, and the same
    // p+p-derived payload is used for both systems; Au+Au inputs are already
    // underlying-event subtracted.
    auto* calibration = new JetCalib("JetCalib_" + suffix.str());
    calibration->set_UseEMfracCalib(false);
    calibration->set_InputNode(collection.rawNode);
    calibration->set_OutputNode(collection.correctedNode);
    calibration->set_JetRadius(static_cast<float>(collection.radius));
    calibration->set_ApplyZvrtxDependentCalib(true);
    calibration->set_ApplyEtaDependentCalib(true);
    server->registerSubsystem(calibration);
  }
}

// Truth jets clustered from generator particles when the input does not ship
// them. No energy scale is ever applied to truth jets.
void registerTruthJets(Fun4AllServer* server, const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  if (!reco.isSimulation || reco.truthJetMode != TruthJetMode::Build)
  {
    return;
  }

  std::set<std::pair<std::string, double>> collections;
  for (const photonjet::JetNodeConfig& collection : plan.tree.jets.nodes)
  {
    if (!collection.truthNode.empty())
    {
      collections.emplace(collection.truthNode, collection.radius);
    }
  }

  for (const auto& [node, radius] : collections)
  {
    auto* truthReco = new JetReco("TruthJetReco_" + node);
    auto* particles = new TruthJetInput(Jet::PARTICLE);
    particles->add_embedding_flag(1);  // the embedded hard-scatter particles only
    truthReco->add_input(particles);
    truthReco->add_algo(antiKt(static_cast<float>(radius), false), node);
    truthReco->set_algo_node("ANTIKT");
    truthReco->set_input_node("TRUTH");
    server->registerSubsystem(truthReco);
  }
}

// Photon candidates. The builder attaches the complete named shower-shape
// family and the isolation sums to every candidate above the storage
// threshold. No identification model is registered here: scoring is a
// downstream stage that reads the stored primitives.
void registerPhotons(Fun4AllServer* server, const Plan& plan, PhotonJetTree* producer)
{
  const Reconstruction& reco = plan.reconstruction;
  if (!reco.reconstructPhotons)
  {
    return;
  }

  // This must run before PhotonClusterBuilder reads its topo-cluster node.
  // absE preserves signed clusters; the photon cone sum uses their signed ET.
  if (reco.useTopoclusterIsolation)
  {
    auto* topo = new RawClusterBuilderTopo("PhotonIsolationTopoClusters");
    topo->set_nodename(reco.topoclusterNode);
    topo->setInputTowerNodePrefix("TOWERINFO_CALIB");
    topo->set_enable_HCal(true);
    topo->set_enable_EMCal(true);
    topo->set_noise(0.0053, 0.0351, 0.0684);
    topo->set_significance(4.0, 2.0, 1.0);
    topo->allow_corner_neighbor(true);
    topo->set_do_split(true);
    topo->set_minE_local_max(1.0, 2.0, 0.5);
    topo->set_R_shower(0.025);
    topo->set_use_only_good_towers(true);
    topo->set_absE(true);
    topo->Verbosity(0);
    server->registerSubsystem(topo);
  }

  auto* photons = new PhotonClusterBuilder("PhotonClusterBuilder");
  photons->set_input_cluster_node(plan.tree.photon.rawClusterNode);
  photons->set_output_photon_node(plan.tree.photon.photonNode);
  photons->set_ET_threshold(static_cast<float>(reco.photonMinEtGeV));
  photons->set_shower_shape_min_tower_energy(static_cast<float>(reco.showerShapeTowerMinEnergyGeV));
  // Shower and isolation thresholds describe different constituent policies.
  // AuAu zero means NO isolation energy cut, including negative SUB1 towers.
  photons->set_isolation_min_tower_energy(reco.isPP ? 0.12F : 0.0F);
  // The standard pp-SIM receipts establish GlobalMbd. Keep the pre-existing
  // archived-DI MBD policy until that separate runtime binding is recovered.
  photons->set_vertex_source(reco.isArchivedDoubleInteraction
      ? PhotonClusterBuilder::VertexSource::Mbd
      : reco.isPP && reco.isSimulation
      ? PhotonClusterBuilder::VertexSource::GlobalMbd
      : PhotonClusterBuilder::VertexSource::MbdThenGlobal);
  photons->set_use_explicit_tower_channel_lookup(true);
  // Both systems use the same configured shower construction (H70 nominal).
  // Additional views use that implementation under named keys. Capture routes
  // by the published floor, so changing a configured view cannot relabel H0/H70.
  for (const std::string& name : plan.tree.showerDefinitions)
  {
    const float floor = name == "H70" ? 0.070F :
                        name == "H0"  ? 0.0F :
                        throw std::runtime_error("unsupported builder shower view: " + name);
    if (floor != static_cast<float>(reco.showerShapeTowerMinEnergyGeV))
    {
      photons->add_shower_shape_view(name, floor);
    }
  }
  photons->set_enable_3x3_moments(true);
  // Preserve the established native COG axis and geometrical support separately
  // from the primary shower-view choice. The AuAu analysis shower is still H70.
  const float nativeAxisFloor = reco.isPP ? 0.070F : 0.0F;
  photons->set_isolation_axis_cog(true, nativeAxisFloor);
  photons->set_require_complete_shower_window(true, nativeAxisFloor);
  photons->set_vertex_cut(reco.photonUseVertexCut, static_cast<float>(reco.photonVertexCutCm));
  photons->set_do_subtracted_iso(reco.useSubtractedIsolation);
  photons->set_do_topocluster_isolation(reco.useTopoclusterIsolation);
  photons->set_topocluster_node(reco.topoclusterNode);
  photons->set_do_bdt(false);
  producer->setPhotonBuilder(photons);
  server->registerSubsystem(photons);
}
}  // namespace

namespace photonjet
{
namespace production
{
using detail::fail;

// RegisterReconstruction
//
void RegisterReconstruction(Fun4AllServer* server, const Plan& plan, PhotonJetTree* producer)
{
  if (!server) fail("null Fun4AllServer");
  if (!producer) fail("null producer");

  const Reconstruction& reco = plan.reconstruction;
  const Conditions& cond = plan.conditions;

  // 0. The exposure observer sees every input event first.
  server->registerSubsystem(new ExposureObserver(producer, plan.tree.nodes.eventHeader));

  // 1. Archived double interaction: replay the frozen random sequence.
  registerArchivedDoubleInteractionSeeds(plan);

  // 2. Minimum-bias detector, zero-degree calorimeter, global vertex.
  registerEventGeometry(server, plan);

  // 3. Calorimeter calibration, status and clusters.
  std::shared_ptr<CemcStatusState> cemcStatus;
  if (cond.cemcStatusGuard)
  {
    cemcStatus = std::make_shared<CemcStatusState>();
    cemcStatus->expectedRun = plan.job.run;
    cemcStatus->originalCalorimeterInput = cond.cemcOriginalCalorimeterInput;
    cemcStatus->recoveryAuthorized = cond.cemcStatusRecoveryAuthorized;
    cemcStatus->recoveryPayload = cond.cemcStatusRecoveryPayload;
  }
  registerCalorimeter(server, plan, cemcStatus);

  // 4. Au+Au data minimum bias and centrality.
  registerCentrality(server, plan);

  // 5. Au+Au underlying-event subtraction.
  registerAuAuBackground(server, plan);

  // 6. Reconstructed jets and the jet energy scale, once each.
  registerRecoJets(server, plan);

  // 7. Truth jets when the input does not ship them.
  registerTruthJets(server, plan);

  // 8. Photon candidates with their shower and isolation primitives.
  registerPhotons(server, plan, producer);

  // 9. The producer, last, so it reads a complete event.
  server->registerSubsystem(producer);

  (void) reco;
}
}  // namespace production
}  // namespace photonjet
