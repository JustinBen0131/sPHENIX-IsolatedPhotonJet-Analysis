/**
 * @file Setup.cc
 * Prepare the detector state consumed by PhotonJetTree.
 * ConfigureConditions binds payloads; RegisterReconstruction builds the event;
 * RegisterInputs opens the exact synchronized source bundle.
 * Configuration parsing stays in Production.cc; photon arithmetic stays in
 * the bundled photonjet::PhotonClusterBuilder.
 */

#include "Production.h"
#include "../EventTree.h"
#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBTTree.h>
#include <calobase/TowerInfoDefs.h>
#include <phool/recoConsts.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>  // Complete type required by MbdReco's unique_ptr<TF1>.
#include <TRandom3.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include "../PhotonJetTree.h"
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
#include <calobase/RawTowerDefs.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <calowaveformsim/CaloWaveformSim.h>
#include "../PhotonClusterBuilder.h"
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterBuilderTopo.h>
#include <calostatusskimmer/CaloStatusSkimmer.h>
#include <calotrigger/MinimumBiasClassifier.h>
#include <centrality/CentralityReco.h>
#include <g4mbd/MbdDigitization.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexReco.h>
#include <mbd/MbdReco.h>
#include <zdcinfo/ZdcReco.h>
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
#include <phool/RunnumberRange.h>
#include <array>
#include <cstdint>
#include <exception>
#include <limits>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>
#include <unordered_map>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/InputFileHandlerReturnCodes.h>

// Conditions and immutable calibration payloads

namespace photonjet
{
namespace production
{
using detail::fail;
using detail::lower;
using detail::trim;
using detail::validSha256;
using detail::fileSha256;
using detail::requireStatus;
using detail::kArchivedDoubleInteractionRun;
using detail::kCemcChannelCount;

namespace
{
bool fileReadable(const std::string& path)
{
  std::ifstream input(path, std::ios::binary);
  return input.good();
}

void verifySha256(const std::string& path, const std::string& expected)
{
  if (path.empty())
  {
    fail("cannot verify an empty payload path");
  }
  if (!fileReadable(path))
  {
    fail("payload is not readable: " + path);
  }
  if (!validSha256(expected))
  {
    fail("invalid expected SHA256 for: " + path);
  }

  const std::string observed = fileSha256(path);
  if (observed != lower(expected))
  {
    fail("SHA256 mismatch for " + path + " expected=" + expected + " observed=" + observed);
  }
}

// The CEMC zero-suppression cross calibration interprets a stored ratio of
// zero as unity and a missing field as NaN. Validate every channel of the
// pinned payload before the calibrator consumes it.
unsigned validateZeroSuppressionPayload(const std::string& path)
{
  requireStatus(!path.empty() && fileReadable(path), "missing zero-suppression payload");
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  requireStatus(file && !file->IsZombie(), "unreadable zero-suppression payload");
  file->Close();

  CDBTTree payload(path);
  unsigned zeroSentinels = 0;
  for (unsigned channel = 0; channel < kCemcChannelCount; ++channel)
  {
    const unsigned key = TowerInfoDefs::encode_emcal(channel);
    const float ratio = payload.GetFloatValue(key, "ratio", 0);
    requireStatus(std::isfinite(ratio) && ratio >= 0,
                  "missing or negative zero-suppression ratio at channel " + std::to_string(channel));
    zeroSentinels += ratio == 0;
  }
  return zeroSentinels;
}

// The run the framework assigns to every input event: the job's own run, the
// data run an embedded simulation was embedded in, or run 28 for archived
// double-interaction hits.
int conditionsRun(const Plan& plan)
{
  if (plan.reconstruction.isArchivedDoubleInteraction) return kArchivedDoubleInteractionRun;
  if (plan.reconstruction.isEmbedded) return static_cast<int>(plan.conditions.cdbTimestamp);
  return plan.job.run;
}
}  // namespace

namespace detail
{
// The file actually resolved at run time is hashed, never a label or a path
// that merely claims to be the payload. Plan loading uses the same hash route.
bool validSha256(const std::string& digest)
{
  return digest.size() == 64 &&
         std::all_of(digest.begin(), digest.end(), [](const unsigned char c) { return std::isxdigit(c) != 0; });
}

std::string fileSha256(const std::string& path) { return eventtree::fileSha256(path); }

void requireStatus(const bool condition, const std::string& reason)
{
  if (!condition)
  {
    throw std::runtime_error("CEMC status contract: " + reason);
  }
}
}  // namespace detail

// ConfigureConditions
//
void ConfigureConditions(Plan& plan)
{
  Conditions& cond = plan.conditions;
  const Reconstruction& reco = plan.reconstruction;

  recoConsts* rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", conditionsRun(plan));
  rc->set_uint64Flag("TIMESTAMP", cond.cdbTimestamp);
  if (!cond.cdbGlobalTag.empty())
  {
    rc->set_StringFlag("CDB_GLOBALTAG", cond.cdbGlobalTag);
  }

  if (reco.reconstructCentrality && cond.centralityLocalPayloads)
  {
    verifySha256(cond.centralityDivisions, cond.centralityDivisionsSha256);
    verifySha256(cond.centralityRunScale, cond.centralityRunScaleSha256);
    verifySha256(cond.centralityVertexScale, cond.centralityVertexScaleSha256);
  }

  if (cond.cemcStatusRecoveryAuthorized)
  {
    verifySha256(cond.cemcStatusRecoveryPayload, cond.cemcStatusRecoveryPayloadSha256);
  }

  if (cond.cemcZeroSuppressionCrossCalibration)
  {
    cond.cemcZeroSuppressionPayload = CDBInterface::instance()->getUrl("CEMC_ZSCrossCalib");
    const unsigned zeroSentinels = validateZeroSuppressionPayload(cond.cemcZeroSuppressionPayload);
    std::cout << "CEMC_ZS_PAYLOAD path=" << cond.cemcZeroSuppressionPayload
              << " zero_means_unity_channels=" << zeroSentinels << std::endl;
  }

  if (cond.applyJetEnergyScale)
  {
    // Resolve the same key the legacy method consumes, and verify those bytes.
    const std::string payload = CDBInterface::instance()->getUrl(cond.jetEnergyScaleCdbKey);
    if (payload.empty())
    {
      fail("the conditions database returned no jet energy-scale payload for key " + cond.jetEnergyScaleCdbKey);
    }
    verifySha256(payload, cond.jetEnergyScalePayloadSha256);
    {
      TDirectory::TContext context;
      std::unique_ptr<TFile> file(TFile::Open(payload.c_str(), "READ"));
      if (!file || file->IsZombie() || !file->IsOpen())
      {
        fail("resolved jet energy-scale payload cannot be opened: " + payload);
      }
    }
    cond.jetEnergyScalePayload = payload;
    plan.tree.jets.energyScalePayloadPath = payload;
  }

  if (cond.applyVertexReweight && !cond.vertexReweightFileSha256.empty())
  {
    verifySha256(cond.vertexReweightFile, cond.vertexReweightFileSha256);
  }

  if (!plan.tree.photon.ppg12TowerMaskFile.empty())
  {
    verifySha256(plan.tree.photon.ppg12TowerMaskFile, plan.tree.photon.ppg12TowerMaskSha256);
  }
}

}  // namespace production
}  // namespace photonjet

// Reconstruction in registration order

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

// The archived double-interaction pedestal events: one of these run-54256
// pedestal files, chosen by the replayed random sequence.
constexpr int kArchivedDoubleInteractionPedestalFiles = 3260;

// Exposure observer
//
// Registered before every other module. It reports each input event to the
// producer before reconstruction has a chance to abort it, which is the only
// way the producer can count the events it never sees. An event whose header
// names another run than the conditions run belongs to a different source and
// stops the job.
//
// ============================================================================
// Source inputs and their read position
//
// A source's DST streams advance together, one source entry per event. Fun4All
// opens a listed file only when the event loop starts, and in this release its
// skip leaves an unopened synchronized input in place and moves a NoSync input
// backwards. A source that starts past its first entry is therefore opened and
// moved here, stream by stream, and every event checks that each stream is at
// the entry the producer records: a misplaced input fails instead of
// relabelling events.
// ============================================================================

class SourceCursor
{
 public:
  virtual ~SourceCursor() = default;

  // Open the stream and move it forward by the given number of entries.
  virtual void advance(std::int64_t entries) = 0;

  // The source entry of the event last read.
  virtual std::int64_t entry() = 0;
};

template <class Manager>
class SourceInput final : public Manager, public SourceCursor
{
 public:
  explicit SourceInput(const std::string& name) : Manager(name) {}

  void advance(const std::int64_t entries) override
  {
    if (!this->IsOpen() && this->OpenNextFile() != InputFileHandlerReturnCodes::SUCCESS)
    {
      throw std::runtime_error("cannot open source input " + this->Name());
    }
    if (!this->IManager() || entries > std::numeric_limits<int>::max() ||
        this->PushBackEvents(-static_cast<int>(entries)) != 0)
    {
      throw std::runtime_error("cannot move source input " + this->Name() + " to its first entry");
    }
  }

  std::int64_t entry() override
  {
    if (!this->IManager())
    {
      throw std::runtime_error("source input " + this->Name() + " is not open");
    }
    return static_cast<std::int64_t>(this->IManager()->getEventNumber()) - 1;
  }
};

std::string inputManagerName(const InputStream& stream) { return "DST_" + stream.name + "_IN"; }

SourceCursor& sourceCursor(const std::string& name)
{
  auto* cursor = dynamic_cast<SourceCursor*>(Fun4AllServer::instance()->getInputManager(name));
  if (!cursor)
  {
    throw std::runtime_error("source input " + name + " is not registered");
  }
  return *cursor;
}

class ExposureObserver final : public SubsysReco
{
 public:
  ExposureObserver(PhotonJetTree* producer, const Plan& plan)
      : SubsysReco("PhotonJetTreeExposureObserver"),
        m_producer(producer),
        m_eventHeaderNode(plan.tree.nodes.eventHeader),
        m_run(photonjet::production::conditionsRun(plan)),
        m_nextEntry(plan.job.firstEntry)
  {
    for (const InputStream& stream : plan.inputs)
    {
      if (!stream.files.empty()) m_inputs.push_back(inputManagerName(stream));
    }
  }

  int process_event(PHCompositeNode* topNode) override
  {
    // Every source stream holds the entry this event is recorded as.
    const std::int64_t expected = m_nextEntry++;
    for (const std::string& name : m_inputs)
    {
      const std::int64_t entry = sourceCursor(name).entry();
      if (entry != expected)
      {
        throw std::runtime_error("source input " + name + " is at entry " + std::to_string(entry) +
                                 " where the event loop expects " + std::to_string(expected));
      }
    }

    std::int64_t sequence = -1;

    if (auto* header = findNode::getClass<EventHeader>(topNode, m_eventHeaderNode); header && header->isValid())
    {
      const int run = header->get_RunNumber();
      if (run > 0 && run != m_run)
      {
        throw std::runtime_error("input event run " + std::to_string(run) + " disagrees with the conditions run " +
                                 std::to_string(m_run));
      }
      sequence = header->get_EvtSequence();
    }

    m_producer->observeUpstreamEvent(sequence);
    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  PhotonJetTree* m_producer;
  std::string m_eventHeaderNode;
  int m_run;
  std::int64_t m_nextEntry;
  std::vector<std::string> m_inputs;
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
// Once FastJetAlgoSub computes areas in coresoftware, this class can be
// retired; until then it stays with the producer that needs it.
//
class JetAreaAlgorithm final : public JetAlgo
{
 public:
  explicit JetAreaAlgorithm(const FastJetOptions& options) : m_options(options)
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
    if (m_options.algo == Jet::ANTIKT)
      os << "ANTIKT";
    else if (m_options.algo == Jet::KT)
      os << "KT";
    else if (m_options.algo == Jet::CAMBRIDGE)
      os << "CAMBRIDGE";
    else
      os << "UNKNOWN";
    os << " r=" << m_options.jet_R << " active_area=" << (m_options.calc_area ? 1 : 0) << std::endl;
  }

  Jet::ALGO get_algo() override { return m_options.algo; }
  float get_par() override { return m_options.jet_R; }

  void cluster_and_fill(std::vector<Jet*>& particles, JetContainer* jetcont) override
  {
    m_fractions.clear();
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
      areaSequence = std::make_unique<fastjet::ClusterSequenceArea>(pseudojets, definition, areaDefinition);
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
      photonjet::JetCalorimeterEnergy calorimeterEnergy;

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
        if (m_options.calc_area && !particle->get_comp_vec().empty())
          calorimeterEnergy.add(static_cast<int>(particle->get_comp_vec().front().first), particle->get_e());
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
      if (m_options.calc_area) m_fractions.emplace(jet, calorimeterEnergy.fractions(totalEnergy));
    }
  }

  std::array<double, 3> fractions(const Jet* jet) const
  {
    const auto found = m_fractions.find(jet);
    if (found == m_fractions.end())
      throw std::runtime_error("raw jet has no current-event calorimeter fraction record");
    return found->second;
  }

 private:
  std::unordered_map<const Jet*, std::array<double, 3>> m_fractions;
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
      return fastjet::JetDefinition(fastjet::kt_algorithm, m_options.jet_R, fastjet::E_scheme, fastjet::Best);
    }
    if (m_options.algo == Jet::CAMBRIDGE)
    {
      return fastjet::JetDefinition(fastjet::cambridge_algorithm, m_options.jet_R, fastjet::E_scheme, fastjet::Best);
    }
    return fastjet::JetDefinition(fastjet::antikt_algorithm, m_options.jet_R, fastjet::E_scheme, fastjet::Best);
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

// JetCalib appends calibrated jets and leaves the output container's default
// metadata untouched. Its output may also be a node loaded from the DST.
// Replace that event's contents and bind its metadata to the raw collection
// before applying the release's calibration, without changing its arithmetic.
class ReconstructedJetCalibration final : public JetCalib
{
 public:
  ReconstructedJetCalibration(const std::string& name, const std::string& inputNode, const std::string& outputNode)
      : JetCalib(name), m_inputNode(inputNode), m_outputNode(outputNode)
  {
    set_InputNode(inputNode);
    set_OutputNode(outputNode);
  }

  int process_event(PHCompositeNode* topNode) override
  {
    auto* raw = findNode::getClass<JetContainer>(topNode, m_inputNode);
    auto* calibrated = findNode::getClass<JetContainer>(topNode, m_outputNode);
    if (!raw || !calibrated || raw == calibrated)
    {
      std::cerr << Name() << ": distinct raw and calibrated jet containers are required" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    calibrated->Reset();
    calibrated->set_algo(raw->get_algo());
    calibrated->set_jetpar_R(raw->get_jetpar_R());
    return JetCalib::process_event(topNode);
  }

 private:
  const std::string m_inputNode;
  const std::string m_outputNode;
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
      : SubsysReco(calibrated ? "CemcStatusBeforeClustering" : "CemcStatusBeforeCalibration"),
        m_state(std::move(state)),
        m_calibrated(calibrated)
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
      // TowerInfoContainer::size() scans the whole TClonesArray; count once.
      const auto towerCount = towers ? towers->size() : 0U;
      requireStatus(towers && towerCount == kCemcChannelCount, "TOWERS_CEMC input missing");

      // Read only the original calorimeter file's RUN metadata: the live RUN
      // node is shared with the paired and geometry inputs and its last
      // loaded conditions URL is not necessarily the calorimeter's.
      requireStatus(!m_state->originalCalorimeterInput.empty(), "missing exact original calorimeter source");
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
        for (unsigned i = 0; i < towerCount; ++i)
        {
          auto* t = towers->get_tower_at_channel(i);
          requireStatus(t && !t->get_isHot(), "missing provenance with existing HOT state requires review");
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
      requireStatus(path.find("/cdb/CEMC_BadTowerMap/") != std::string::npos && path.size() > suffix.size() &&
                        path.compare(path.size() - suffix.size(), suffix.size(), suffix) == 0,
                    "payload path is not bound to this run and the CEMC domain");

      std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
      requireStatus(file && !file->IsZombie(), "unreadable CEMC map file");
      file->Close();

      CDBTTree map(path);
      m_state->mask.clear();
      m_state->mask.reserve(kCemcChannelCount);
      for (unsigned i = 0; i < towerCount; ++i)
      {
        const int value = map.GetIntValue(towers->encode_key(i), "status", 0);
        // A missing integer field reads as INT_MIN and must never count as good.
        requireStatus(value >= 0 && value <= 255, "missing or invalid map status field");
        m_state->mask.push_back(value > 0);
      }
      m_state->initialized = true;

      std::cout << "CEMC_STATUS_CONTRACT_INIT run=" << m_state->expectedRun
                << " original_input=" << m_state->originalCalorimeterInput << " payload=" << path
                << " action=" << (m_state->restoringMissingStatus ? "RESTORE_MISSING_HOT_ONLY" : "PRESERVE_VALIDATE")
                << " map_rejected=" << std::count(m_state->mask.begin(), m_state->mask.end(), 1) << std::endl;
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
      auto* towers = findNode::getClass<TowerInfoContainer>(top, m_calibrated ? "TOWERINFO_CALIB_CEMC" : "TOWERS_CEMC");
      // TowerInfoContainer::size() scans the whole TClonesArray; count once
      // per event. Re-counting in these loop conditions doubled AuAu job time
      // in the reference production.
      const auto towerCount = towers ? towers->size() : 0U;
      requireStatus(towers && towerCount == kCemcChannelCount, "required CEMC node missing");

      if (!m_calibrated && m_state->restoringMissingStatus)
      {
        for (unsigned i = 0; i < towerCount; ++i)
        {
          if (m_state->mask[i])
          {
            towers->get_tower_at_channel(i)->set_isHot(true);
          }
        }
      }

      // Same checks, same order and messages; the message strings are built
      // only on failure instead of for all 24,576 channels every event.
      for (unsigned i = 0; i < towerCount; ++i)
      {
        auto* t = towers->get_tower_at_channel(i);
        if (t == nullptr)
        {
          requireStatus(false, "null tower");
        }
        if (m_state->mask[i] && !(t->get_status() & kCemcHotBit))
        {
          requireStatus(false, "map-rejected tower lacks HOT flag at channel " + std::to_string(i));
        }
        if (m_state->mask[i] && !(t->get_isHot() && !t->get_isGood()))
        {
          requireStatus(false, "map-rejected tower must be HOT and not GOOD in the actual node");
        }
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
      std::cout << "CEMC_STATUS_CONTRACT_END run=" << m_state->expectedRun << " raw_events=" << m_state->rawEvents
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

  if (reco.reconstructZdc)
  {
    auto* zdc = new ZdcReco();
    zdc->set_zdc1_cut(0.0);
    zdc->set_zdc2_cut(0.0);
    server->registerSubsystem(zdc);
  }

  server->registerSubsystem(new GlobalVertexReco());
}

// The archived double-interaction reconstruction replays its source's own
// five-seed random sequence, so the reconstructed calorimeter matches the
// archived production exactly. The seeds are consumed first in, first out:
// the MBD digitisation (constructed with the event geometry), the pedestal
// draw, then the three waveform simulations at InitRun.
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
  if (reco.historicalSeeds.size() != 5)
  {
    fail("archived double-interaction reconstruction requires its source's five-seed sequence");
  }

  PHRandomSeed::Verbosity(1);
  for (const unsigned int seed : reco.historicalSeeds)
  {
    PHRandomSeed::LoadSeed(seed);
  }
}

// The archived double-interaction calorimeter: the run-28 standard tower
// helpers (CEMC_Towers, HCALInner_Towers and HCALOuter_Towers of the release
// macros, tower information on, pp zero-angle beam) turn the archived Geant4
// hits into TOWERS_* waveform towers over pedestal events from the file the
// replayed sequence chooses. The ordinary calibration chain follows them.
void registerArchivedDoubleInteractionTowers(Fun4AllServer* server, const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  if (!reco.isArchivedDoubleInteraction)
  {
    return;
  }

  TRandom3 draw;
  draw.SetSeed(PHRandomSeed());
  const int sequence = static_cast<int>(draw.Integer(kArchivedDoubleInteractionPedestalFiles));
  if (sequence != reco.expectedPedestalSequence)
  {
    fail("the replayed seeds chose pedestal file " + std::to_string(sequence) + ", not the archived " +
         std::to_string(reco.expectedPedestalSequence));
  }
  std::ostringstream pedestalFile;
  pedestalFile << "pedestal-54256-0" << std::setw(4) << std::setfill('0') << sequence << ".root";
  auto* pedestals = new Fun4AllNoSyncDstInputManager("DST2");
  pedestals->AddFile(pedestalFile.str());
  pedestals->Repeat();
  server->registerInputManager(pedestals);

  const char* calibrationRoot = std::getenv("CALIBRATIONROOT");
  if (!calibrationRoot || !*calibrationRoot)
  {
    fail("CALIBRATIONROOT is required for the CEMC light-collection model");
  }

  struct Calorimeter
  {
    std::string name;
    CaloTowerDefs::DetectorSystem type;
    int zeroSuppressionAdc;
  };
  const Calorimeter calorimeters[] = {{"CEMC", CaloTowerDefs::CEMC, 60},
                                      {"HCALIN", CaloTowerDefs::HCALIN, 30},
                                      {"HCALOUT", CaloTowerDefs::HCALOUT, 30}};
  for (const Calorimeter& calorimeter : calorimeters)
  {
    auto* waveforms = new CaloWaveformSim(calorimeter.name + "CaloWaveformSim");
    waveforms->set_detector_type(calorimeter.type);
    waveforms->set_detector(calorimeter.name);
    waveforms->set_timewidth(0.2);
    waveforms->set_peakpos(6);
    if (calorimeter.type == CaloTowerDefs::CEMC)
    {
      waveforms->get_light_collection_model().load_data_file(
          std::string(calibrationRoot) + "/CEMC/LightCollection/Prototype3Module.xml",
          "data_grid_light_guide_efficiency", "data_grid_fiber_trans");
      waveforms->set_pedestal_scale(0.75);
    }
    if (calorimeter.type == CaloTowerDefs::HCALOUT)
    {
      waveforms->set_calibName("HCALOUT_calib_ADC_to_ETower");
    }
    server->registerSubsystem(waveforms);

    auto* towers = new CaloTowerBuilder(calorimeter.name + "CaloTowerBuilder");
    towers->set_detector_type(calorimeter.type);
    towers->set_nsamples(12);
    towers->set_dataflag(false);
    towers->set_processing_type(CaloWaveformProcessing::TEMPLATE);
    towers->set_builder_type(CaloTowerDefs::kWaveformTowerSimv1);
    towers->set_softwarezerosuppression(true, calorimeter.zeroSuppressionAdc);
    server->registerSubsystem(towers);

    auto* status = new CaloTowerStatus(calorimeter.name + "STATUS");
    status->set_detector_type(calorimeter.type);
    server->registerSubsystem(status);

    auto* calibration = new CaloTowerCalib(calorimeter.name + "CALIB");
    calibration->set_detector_type(calorimeter.type);
    calibration->set_outputNodePrefix("TOWERINFO_CALIB_");
    server->registerSubsystem(calibration);
  }
}

// Calorimeter calibration, status and clustering, in the order of the
// historical Process_Calo_Calib.
void registerCalorimeter(Fun4AllServer* server, const Plan& plan, const std::shared_ptr<CemcStatusState>& cemcStatus)
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
    const bool fastProcessing = (timestamp > RunnumberRange::RUN2PP_FIRST && timestamp < RunnumberRange::RUN2PP_LAST) ||
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

  // Raw-tower status validation against the original bad-tower map.
  if (cemcStatus)
  {
    server->registerSubsystem(new CemcTowerStatusGuard(cemcStatus, false));
  }

  // Tower status setters: simulation only (see Production.cc).
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
void registerRecoJets(Fun4AllServer* server, const Plan& plan, PhotonJetTree* producer)
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
    auto* algorithm = antiKt(static_cast<float>(collection.radius), true);
    jetReco->add_algo(algorithm, collection.rawNode);
    producer->setJetFractionReader(collection.rawNode,
                                   [algorithm](const Jet* jet) { return algorithm->fractions(jet); });
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
    auto* calibration =
        new ReconstructedJetCalibration("JetCalib_" + suffix.str(), collection.rawNode, collection.correctedNode);
    calibration->set_UseEMfracCalib(false);
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
// threshold. The configured named models (the PPG12 NPB and producer scores)
// only add scores, each bound to its file digest; no model selects or removes
// a candidate, and the nominal photon ID is filled downstream.
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

  auto* photons = new photonjet::PhotonClusterBuilder("PhotonClusterBuilder");
  photons->set_capture_model_inputs(true);
  photons->set_input_cluster_node(plan.tree.photon.rawClusterNode);
  photons->set_output_photon_node(plan.tree.photon.photonNode);
  photons->set_ET_threshold(static_cast<float>(reco.photonMinEtGeV));
  photons->set_shower_shape_min_tower_energy(static_cast<float>(reco.showerShapeTowerMinEnergyGeV));
  // Shower and isolation thresholds describe different constituent policies.
  // AuAu zero means NO isolation energy cut, including negative SUB1 towers.
  photons->set_isolation_min_tower_energy(reco.isPP ? 0.12F : 0.0F);
  // The standard pp-SIM receipts establish GlobalMbd. Keep the pre-existing
  // archived-DI MBD policy until that separate runtime binding is recovered.
  photons->set_vertex_source(reco.isArchivedDoubleInteraction ? photonjet::PhotonClusterBuilder::VertexSource::Mbd
                             : reco.isPP && reco.isSimulation
                                 ? photonjet::PhotonClusterBuilder::VertexSource::GlobalMbd
                                 : photonjet::PhotonClusterBuilder::VertexSource::MbdThenGlobal);
  photons->set_use_explicit_tower_channel_lookup(true);
  // Seven frozen views use one configurable builder algorithm. Primary H70
  // outputs, candidate support and isolation remain independent of audit views.
  for (const std::string& name : plan.tree.showerDefinitions)
  {
    // The view table is the single definition of each retained view; its
    // semantic digest is recorded with the output.
    const photonjet::ShowerViewDefinition& view = photonjet::showerViewDefinition(name);
    const float floor = static_cast<float>(view.floorGeV);
    if (name == "H0" && floor != static_cast<float>(reco.showerShapeTowerMinEnergyGeV))
    {
      photons->add_shower_shape_view(name, floor);
    }
    using Builder = photonjet::PhotonClusterBuilder;
    const auto membership = [](std::int32_t owned)
    { return owned ? Builder::ShowerMembership::Owned : Builder::ShowerMembership::Grid; };
    photons->add_shower_feature_view(
        name, floor, view.energySource ? Builder::ShowerEnergy::ClusterMap : Builder::ShowerEnergy::Calibrated,
        membership(view.sumMembership), membership(view.momentMembership));
  }
  // The established Au+Au productions also gave each photon its 3x3 core
  // moments (pho_weta33_cogx, pho_wphi33_cogx, model inputs); p+p carried them
  // only through the shower views.
  photons->set_enable_3x3_moments(reco.isAuAu);
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
  for (const auto& model : plan.tree.photon.models)
  {
    // Bind actual file bytes before constructing the evaluator. Configured
    // filenames alone are never treated as model provenance.
    if (model.sha256.empty() || photonjet::production::detail::fileSha256(model.file) != model.sha256)
      throw std::runtime_error("photon model file digest differs: " + model.scoreName);
    photons->add_bdt_model(model.scoreName, model.file, model.features);
    if (!std::isnan(model.minEt) || !std::isnan(model.maxEt) || !std::isnan(model.maxAbsEta))
      photons->set_bdt_model_domain(model.scoreName, model.minEt, model.maxEt, model.maxAbsEta);
    photons->set_bdt_model_unevaluated_score(model.scoreName, model.unevaluatedScore);
  }
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
  server->registerSubsystem(new ExposureObserver(producer, plan));

  // 1. Archived double interaction: replay the frozen random sequence.
  registerArchivedDoubleInteractionSeeds(plan);

  // 2. Archived hits need their reconstructed vertex before calorimetry, and
  // their towers built from the hits over replayed pedestal events. Ordinary
  // DATA clustering instead uses the vertex shipped on the DST, exactly as
  // Process_Calo_Calib did before the later MBD reconstruction.
  if (reco.isArchivedDoubleInteraction)
  {
    registerEventGeometry(server, plan);
    registerArchivedDoubleInteractionTowers(server, plan);
  }

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

  // Preserve the original registration order: clusters first, then MBD,
  // ZDC and GlobalVertex, then photons/jets using that refreshed vertex.
  // Moving MBD ahead of clustering changes the cluster position correction.
  if (!reco.isArchivedDoubleInteraction) registerEventGeometry(server, plan);

  // 4. Au+Au data minimum bias and centrality.
  registerCentrality(server, plan);

  // 5. Au+Au underlying-event subtraction.
  registerAuAuBackground(server, plan);

  // 6. Reconstructed jets and the jet energy scale, once each.
  registerRecoJets(server, plan, producer);

  // 7. Truth jets when the input does not ship them.
  registerTruthJets(server, plan);

  // 8. Photon candidates with their shower and isolation primitives.
  registerPhotons(server, plan, producer);

  // 9. The producer, last, so it reads a complete event.
  server->registerSubsystem(producer);
}
}  // namespace production
}  // namespace photonjet

// Physical streams and synchronized input registration

namespace photonjet
{
namespace production
{

using detail::fail;
using detail::trim;

namespace
{

// ============================================================================
// Input-list row parsing
//
// The frontend binds and hashes the input-list file. This layer alone
// interprets its physical column layout.
//
// Empty-lane sentinels are normalized to empty strings so optional lanes can be
// handled explicitly rather than through filename heuristics.
// ============================================================================

std::vector<std::string> splitColumns(const std::string& line)
{
  std::istringstream stream(line);

  std::vector<std::string> columns;
  std::string token;

  while (stream >> token)
  {
    if (token == "NONE" || token == "none" || token == "-")
    {
      columns.emplace_back();
    }
    else
    {
      columns.push_back(token);
    }
  }

  return columns;
}


/**
 * Append one positional input-list column when that lane is present.
 *
 * Missing optional lanes remain absent rather than being represented by a
 * placeholder filename.
 */
void appendColumn(std::vector<std::string>& destination, const std::vector<std::string>& columns,
                  const std::size_t index)
{
  if (index < columns.size() && !columns[index].empty())
  {
    destination.push_back(columns[index]);
  }
}


// ============================================================================
// Fun4All input-manager construction
//
// InputStream::synchronized is already a resolved production decision. This
// function translates that decision directly into the corresponding Fun4All
// manager type without reinterpreting the source mode.
// ============================================================================

Fun4AllInputManager* makeInputManager(const InputStream& stream, const bool recomputeCentrality)
{
  const std::string name = inputManagerName(stream);

  if (stream.synchronized)
  {
    auto* manager = new SourceInput<Fun4AllDstInputManager>(name);
    // Au+Au DATA recomputes centrality: the paired DSTs' own CentralityInfo
    // is not imported, so CentralityReco creates the node this pass reads.
    // Their MBD inputs are kept.
    if (recomputeCentrality)
    {
      manager->BranchSelect("*CentralityInfo", 0);
    }
    return manager;
  }

  return new SourceInput<Fun4AllNoSyncDstInputManager>(name);
}

}  // namespace


namespace detail
{

// ============================================================================
// Physical input-list interpretation
//
// Convert one source-bound input-list row into named InputStream objects.
//
// This is intentionally the only code that knows:
//
//   - how many columns each InputLayout requires
//   - which physical column corresponds to which logical stream
//   - which lanes are required versus optional
//   - which modes use ordinary Fun4All synchronization
//
// The resulting InputStream inventory is subsequently consumed by
// Production.cc and RegisterInputs without any positional-column knowledge.
// ============================================================================

std::vector<InputStream> loadInputs(const Job& job, const Reconstruction& reconstruction)
{
  std::ifstream input(job.inputList);

  if (!input)
  {
    fail("cannot read input list: " + job.inputList);
  }

  std::vector<std::string> primary;
  std::vector<std::string> auxiliary;
  std::vector<std::string> g4;
  std::vector<std::string> truthJets;
  std::vector<std::string> global;
  std::vector<std::string> mbd;

  std::string line;
  std::size_t sourceRows = 0;

  while (std::getline(input, line))
  {
    line = trim(line);

    if (line.empty() || line.front() == '#')
    {
      continue;
    }

    const std::vector<std::string> columns = splitColumns(line);

    if (columns.empty())
    {
      continue;
    }

    /*
     * One invocation owns one physical source bundle and ultimately one
     * SourceRecord. Multiple source rows would require explicit source-boundary
     * accounting before they could safely share the same producer invocation.
     */
    if (++sourceRows > 1)
    {
      fail("one physical source bundle per invocation is required");
    }

    const std::size_t expected = reconstruction.inputLayout == InputLayout::ReconstructedSimulation ? 5 : 2;

    if (columns.size() != expected)
    {
      fail("input-list row has an unexpected column count");
    }

    // ------------------------------------------------------------------------
    // Positional row -> named logical lanes
    // ------------------------------------------------------------------------

    switch (reconstruction.inputLayout)
    {
      case InputLayout::PairedData:
        appendColumn(primary, columns, 0);

        appendColumn(auxiliary, columns, 1);

        break;

      case InputLayout::ReconstructedSimulation:
        appendColumn(primary, columns, 0);

        appendColumn(g4, columns, 1);

        appendColumn(truthJets, columns, 2);

        appendColumn(global, columns, 3);

        appendColumn(mbd, columns, 4);

        break;

      case InputLayout::ArchivedDoubleInteraction:
        appendColumn(g4, columns, 0);

        appendColumn(truthJets, columns, 1);

        break;
    }
  }

  // ==========================================================================
  // Build the resolved stream inventory
  //
  // Required streams are retained even if their file vector is empty so later
  // Plan validation can fail with an explicit missing-stream error. Optional
  // streams that are absent are omitted entirely.
  // ==========================================================================

  std::vector<InputStream> streams;

  auto add = [&](const char* name, std::vector<std::string> files, const bool requiredStream, const bool synchronized)
  {
    if (files.empty() && !requiredStream)
    {
      return;
    }

    InputStream stream;

    stream.name = name;

    stream.files = std::move(files);

    stream.required = requiredStream;

    stream.synchronized = synchronized;

    streams.push_back(std::move(stream));
  };

  /*
   * Embedded and archived inputs intentionally retain the historical NoSync
   * behavior. Paired DATA and ordinary reconstructed simulation use normal
   * Fun4All synchronization.
   */
  const bool sync = !(reconstruction.isEmbedded || reconstruction.isArchivedDoubleInteraction);

  switch (reconstruction.inputLayout)
  {
    case InputLayout::PairedData:
      add("PRIMARY", std::move(primary), true, true);

      add("AUXILIARY", std::move(auxiliary), false, true);

      break;

    case InputLayout::ArchivedDoubleInteraction:
      add("G4HITS", std::move(g4), true, false);

      add("TRUTH_JETS", std::move(truthJets), true, false);

      break;

    case InputLayout::ReconstructedSimulation:
      add("CALO", std::move(primary), true, sync);

      add("G4HITS", std::move(g4), true, sync);

      add("TRUTH_JETS", std::move(truthJets), false, sync);

      add("GLOBAL", std::move(global), false, sync);

      add("MBD", std::move(mbd), false, sync);

      break;
  }

  return streams;
}


/**
 * Return whether a named resolved input stream exists and carries at least one
 * physical file.
 *
 * Production.cc uses this after list interpretation to resolve modes such as
 * TruthJetMode::Auto without inspecting the positional input-list layout.
 */
bool hasInputStream(const Plan& plan, const std::string& name)
{
  return std::any_of(plan.inputs.begin(), plan.inputs.end(),
                     [&](const InputStream& stream) { return stream.name == name && !stream.files.empty(); });
}

}  // namespace detail


// ============================================================================
// Fun4All input registration
//
// RegisterInputs performs the final translation from named InputStream objects
// into Fun4All input managers.
//
// If calorimeter reconstruction occurs during this pass, the calorimeter
// geometry run node is registered first. Physical source streams are then
// attached in their already-resolved order and synchronization mode.
//
// Detector reconstruction itself is registered separately by
// RegisterReconstruction.
// ============================================================================

void RegisterInputs(Fun4AllServer* server, const Plan& plan)
{
  if (!server)
  {
    fail("null Fun4AllServer");
  }

  // ==========================================================================
  // Calorimeter geometry run node
  //
  // Reconstructing calibrated towers/clusters during this pass requires the
  // calorimeter geometry payload to be available as a run-node input.
  // ==========================================================================

  if (plan.reconstruction.reconstructCalorimeter)
  {
    auto* geometry = new Fun4AllRunNodeInputManager("DST_GEO");

    geometry->AddFile(CDBInterface::instance()->getUrl("calo_geo"));

    server->registerInputManager(geometry);
  }

  // ==========================================================================
  // Physical source streams
  // ==========================================================================

  for (const InputStream& stream : plan.inputs)
  {
    if (stream.files.empty())
    {
      if (stream.required)
      {
        fail("required input stream has no files: " + stream.name);
      }

      continue;
    }

    Fun4AllInputManager* manager = makeInputManager(stream, plan.reconstruction.reconstructCentrality);

    manager->Verbosity(plan.job.verbosity);

    for (const std::string& file : stream.files)
    {
      manager->AddFile(file);
    }

    server->registerInputManager(manager);
  }
}


void StartAtFirstEntry(Fun4AllServer* server, const Plan& plan)
{
  if (!server)
  {
    fail("StartAtFirstEntry requires a Fun4AllServer");
  }

  if (plan.job.firstEntry == 0)
  {
    return;
  }

  for (const InputStream& stream : plan.inputs)
  {
    if (!stream.files.empty())
    {
      sourceCursor(inputManagerName(stream)).advance(plan.job.firstEntry);
    }
  }
}

}  // namespace production
}  // namespace photonjet
