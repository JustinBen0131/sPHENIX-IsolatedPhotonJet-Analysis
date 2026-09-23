//
// Production.cc
//
// Builds the reconstructed event that PhotonJetTree reads.
//
// Owns
//   configuration parsing and profile resolution
//   physical input-list interpretation and Fun4All input managers
//   run context, conditions-database payloads and their digest checks
//   the exposure observer that accounts for events reconstruction rejects
//   calorimeter calibration, tower status and clustering
//   minimum-bias detector, global vertex, zero-degree calorimeter
//   Au+Au minimum-bias classification and centrality
//   Au+Au underlying-event subtraction
//   jet clustering with active area, and the jet energy scale, applied once
//   truth-jet clustering when the input does not ship truth jets
//   photon candidate reconstruction, with no identification model attached
//
// Does not own
//   any event observable, record or tree
//
// Every module below is registered in the order the historical production
// registered it. Where a parameter came from an environment variable in that
// production it is now a typed field of the plan, with the same default.
//
#include "Production.h"

#include "../PhotonJetTree.h"

// ---- sPHENIX framework -----------------------------------------------------

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/CdbUrlSave.h>
#include <ffaobjects/EventHeader.h>
#include <ffaobjects/RunHeader.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
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

// ---- YAML, ROOT, system -------------------------------------------------------

#include <yaml-cpp/yaml.h>

#include <TDirectory.h>
#include <TFile.h>

#include <phool/RunnumberRange.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
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

//
// Frozen production constants
//
// These are the values the reference production accepted. They are checked
// against the configuration rather than read from it, so a configuration
// edit cannot silently move a nominal choice.
//
constexpr const char* kNominalCentralitySet = "sam_fit_20260913_v1";

constexpr const char* kNominalCentralityManifestSha256 =
    "4a11e0e10ae9b82062325c3275679935da2259f4c91f796950a6b5ddd2e1501c";

constexpr const char* kNominalJetEnergyScaleSha256 =
    "76a8788fdb4e0b3859d01361f884e295ea609bb8c1cc188105389214db2de27d";

constexpr int kQualifiedCentralityFallbackRun = 68144;

constexpr int kArchivedDoubleInteractionRun = 28;
constexpr int kArchivedDoubleInteractionPedestalSequence = 534;

constexpr std::array<unsigned int, 5> kArchivedDoubleInteractionSeeds =
{
  2991264730U, 4256268992U, 2394322166U, 874466025U, 2240380304U
};

// Runs at or below this number are simulation timestamps.
constexpr int kDataRunThreshold = 1000;

//
// Small helpers
//
[[noreturn]] void fail(const std::string& message)
{
  throw std::runtime_error("PhotonJetTree production: " + message);
}

std::string trim(const std::string& value)
{
  const auto first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
  {
    return {};
  }
  const auto last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

std::string lower(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(),
                 [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

std::string readTextFile(const std::string& path)
{
  std::ifstream input(path, std::ios::binary);
  if (!input)
  {
    fail("cannot read file: " + path);
  }
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

bool fileReadable(const std::string& path)
{
  std::ifstream input(path, std::ios::binary);
  return input.good();
}

// ${NAME} expansion for configured paths. An unset variable is an error, not
// an empty string, because an empty path would silently select nothing.
std::string expandEnvironment(const std::string& input)
{
  std::string output;
  for (std::size_t i = 0; i < input.size();)
  {
    if (input[i] == '$' && i + 1 < input.size() && input[i + 1] == '{')
    {
      const std::size_t close = input.find('}', i + 2);
      if (close == std::string::npos)
      {
        fail("unterminated environment variable in path: " + input);
      }
      const std::string key = input.substr(i + 2, close - i - 2);
      const char* value = std::getenv(key.c_str());
      if (!value)
      {
        fail("required environment variable is unset: " + key);
      }
      output += value;
      i = close + 1;
      continue;
    }
    output += input[i];
    ++i;
  }
  return output;
}

//
// SHA256 verification
//
// The file actually resolved at run time is hashed, never a label or a path
// that merely claims to be the payload.
//
bool validSha256(const std::string& digest)
{
  return digest.size() == 64 &&
         std::all_of(digest.begin(), digest.end(),
                     [](const unsigned char c) { return std::isxdigit(c) != 0; });
}

bool safeHashPath(const std::string& path)
{
  return !path.empty() &&
         std::all_of(path.begin(), path.end(), [](const unsigned char c) {
           return std::isalnum(c) || c == '/' || c == '.' || c == '_' ||
                  c == '-' || c == '+';
         });
}

std::string fileSha256(const std::string& path)
{
  if (!safeHashPath(path))
  {
    fail("unsafe path supplied for SHA256 verification: " + path);
  }

#ifdef __APPLE__
  const std::string command = "/usr/bin/shasum -a 256 -- " + path;
#else
  const std::string command = "sha256sum -- " + path;
#endif

  FILE* pipe = ::popen(command.c_str(), "r");
  if (!pipe)
  {
    fail("failed to launch the SHA256 tool for: " + path);
  }

  char buffer[256] = {};
  const bool readOk = std::fgets(buffer, sizeof(buffer), pipe) != nullptr;

  errno = 0;
  const int status = ::pclose(pipe);
  const int closeErrno = errno;

  if (!readOk)
  {
    fail("failed to read the SHA256 result for: " + path);
  }

  std::string digest(buffer);
  digest = lower(trim(digest.substr(0, digest.find_first_of(" \t\r\n"))));

  if (!validSha256(digest))
  {
    fail("invalid SHA256 result for: " + path);
  }

  // ROOT or Fun4All can reap the short-lived child before pclose waits for
  // it. That is not a hash failure; the digest comparison still decides.
  const bool externallyReaped = status == -1 && closeErrno == ECHILD;
  if (status != 0 && !externallyReaped)
  {
    fail("SHA256 tool failed for: " + path);
  }

  return digest;
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
    fail("SHA256 mismatch for " + path + " expected=" + expected +
         " observed=" + observed);
  }
}

//
// YAML helpers
//
template <class T>
T required(const YAML::Node& node, const std::string& key, const std::string& context)
{
  if (!node || !node[key] || node[key].IsNull())
  {
    fail("missing required configuration key '" + key + "' in " + context);
  }
  return node[key].as<T>();
}

template <class T>
T optional(const YAML::Node& node, const std::string& key, const T& fallback)
{
  if (!node || !node[key] || node[key].IsNull())
  {
    return fallback;
  }
  return node[key].as<T>();
}

std::string optionalString(const YAML::Node& node, const std::string& key,
                           const std::string& fallback = {})
{
  if (!node || !node[key] || node[key].IsNull())
  {
    return fallback;
  }
  return node[key].as<std::string>();
}

std::vector<std::string> stringList(const YAML::Node& node)
{
  std::vector<std::string> result;
  if (node && node.IsSequence())
  {
    for (const YAML::Node& value : node)
    {
      result.push_back(value.as<std::string>());
    }
  }
  return result;
}

std::vector<double> doubleList(const YAML::Node& node)
{
  std::vector<double> result;
  if (node && node.IsSequence())
  {
    for (const YAML::Node& value : node)
    {
      result.push_back(value.as<double>());
    }
  }
  return result;
}

std::vector<int> intList(const YAML::Node& node)
{
  std::vector<int> result;
  if (node && node.IsSequence())
  {
    for (const YAML::Node& value : node)
    {
      result.push_back(value.as<int>());
    }
  }
  return result;
}

TruthJetMode parseTruthJetMode(const std::string& value)
{
  const std::string mode = lower(value);
  if (mode.empty() || mode == "none") return TruthJetMode::None;
  if (mode == "auto") return TruthJetMode::Auto;
  if (mode == "dst") return TruthJetMode::DST;
  if (mode == "build") return TruthJetMode::Build;
  fail("unknown truth-jet mode: " + value);
}

photonjet::JetView parseJetView(const std::string& value)
{
  const std::string view = lower(value);
  if (view == "pp") return photonjet::JetView::PP;
  if (view == "auau_sub1" || view == "sub1") return photonjet::JetView::AuAuSub1;
  if (view == "auau_nosub" || view == "nosub") return photonjet::JetView::AuAuNoSub;
  fail("unknown jet view: " + value);
}

//
// Input-list parsing
//
// One line per source. Columns are whitespace separated and "NONE" or "-"
// marks an absent lane. This is the only code that knows the column layout.
//
//   paired data                     CALO  AUXILIARY
//   reconstructed simulation        CALO  G4  TRUTH_JETS  GLOBAL  MBD
//   archived double interaction     G4  TRUTH_JETS
//
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

void appendColumn(std::vector<std::string>& destination,
                  const std::vector<std::string>& columns,
                  const std::size_t index)
{
  if (index < columns.size() && !columns[index].empty())
  {
    destination.push_back(columns[index]);
  }
}

std::vector<InputStream> loadInputs(const Job& job, const Reconstruction& reconstruction)
{
  std::ifstream input(job.inputList);
  if (!input)
  {
    fail("cannot read input list: " + job.inputList);
  }

  std::vector<std::string> primary, auxiliary, g4, truthJets, global, mbd;

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

    if (++sourceRows > 1) fail("one physical source bundle per invocation is required");
    const std::size_t expected = reconstruction.inputLayout == InputLayout::ReconstructedSimulation ? 5 : 2;
    if (columns.size() != expected) fail("input-list row has an unexpected column count");
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

  std::vector<InputStream> streams;

  auto add = [&](const char* name, std::vector<std::string> files,
                 const bool requiredStream, const bool synchronized)
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

  // Embedded and archived inputs are read without Fun4All synchronisation,
  // as the historical production did; paired data is synchronised.
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

bool hasInputStream(const Plan& plan, const std::string& name)
{
  return std::any_of(plan.inputs.begin(), plan.inputs.end(),
                     [&](const InputStream& stream) {
                       return stream.name == name && !stream.files.empty();
                     });
}

//
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
constexpr std::size_t kCemcChannelCount = 24576;
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

void requireStatus(const bool condition, const std::string& reason)
{
  if (!condition)
  {
    throw std::runtime_error("CEMC status contract: " + reason);
  }
}

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
                  "missing or negative zero-suppression ratio at channel " +
                      std::to_string(channel));
    zeroSentinels += ratio == 0;
  }
  return zeroSentinels;
}

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
void registerPhotons(Fun4AllServer* server, const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  if (!reco.reconstructPhotons)
  {
    return;
  }

  auto* photons = new PhotonClusterBuilder("PhotonClusterBuilder");
  photons->set_input_cluster_node(plan.tree.photon.rawClusterNode);
  photons->set_output_photon_node(plan.tree.photon.photonNode);
  photons->set_ET_threshold(static_cast<float>(reco.photonMinEtGeV));
  photons->set_shower_shape_min_tower_energy(static_cast<float>(reco.showerShapeTowerMinEnergyGeV));
  photons->set_enable_3x3_moments(true);
  photons->set_isolation_axis_cog(true);
  photons->set_vertex_cut(reco.photonUseVertexCut, static_cast<float>(reco.photonVertexCutCm));
  photons->set_do_subtracted_iso(reco.useSubtractedIsolation);
  photons->set_do_topocluster_isolation(reco.useTopoclusterIsolation);
  photons->set_topocluster_node(reco.topoclusterNode);
  photons->set_do_bdt(false);
  server->registerSubsystem(photons);
}

Fun4AllInputManager* makeInputManager(const InputStream& stream)
{
  const std::string name = "DST_" + stream.name + "_IN";
  if (stream.synchronized)
  {
    return new Fun4AllDstInputManager(name);
  }
  return new Fun4AllNoSyncDstInputManager(name);
}

}  // namespace

//
// Public interface
//
namespace photonjet
{
namespace production
{

//
// LoadPlan
//
Plan LoadPlan(const std::string& yamlFile, const Job& job)
{
  if (job.profile.empty()) fail("production profile is empty");
  if (job.inputList.empty()) fail("input-list path is empty");
  if (job.outputFile.empty()) fail("output file name is empty");

  Plan plan;
  plan.job = job;

  plan.configurationText = readTextFile(yamlFile);
  plan.configurationSha256 = fileSha256(yamlFile);

  const YAML::Node root = YAML::Load(plan.configurationText);
  const YAML::Node profile = root["profiles"][job.profile];
  if (!profile)
  {
    fail("profile '" + job.profile + "' is absent from " + yamlFile);
  }
  const std::string context = "profiles." + job.profile;

  Reconstruction& reco = plan.reconstruction;
  Conditions& cond = plan.conditions;
  photonjet::Config& tree = plan.tree;

  // ---- mode ---------------------------------------------------------------

  const std::string system = lower(required<std::string>(profile, "system", context));
  const std::string dataKind = lower(required<std::string>(profile, "data_kind", context));
  const std::string role = lower(required<std::string>(profile, "simulation_role", context));
  const std::string inputMode = lower(required<std::string>(profile, "input_mode", context));

  reco.isPP = system == "pp";
  reco.isAuAu = system == "auau";
  reco.isData = dataKind == "data";
  reco.isSimulation = dataKind == "simulation";
  reco.isEmbedded = reco.isAuAu && reco.isSimulation;
  reco.isArchivedDoubleInteraction = inputMode == "archived_g4_only";

  if (reco.isPP == reco.isAuAu) fail("profile must name exactly one collision system");
  if (reco.isData == reco.isSimulation) fail("profile must name exactly one data kind");

  tree.system = reco.isPP ? CollisionSystem::PP : CollisionSystem::AuAu;
  tree.dataKind = reco.isData ? DataKind::Data : DataKind::Simulation;
  tree.inputMode = reco.isArchivedDoubleInteraction ? InputMode::ArchivedG4Only
                                                    : InputMode::ReconstructedDST;
  if (role == "photon_jet") tree.simulationRole = SimulationRole::PhotonJet;
  else if (role == "inclusive_jet") tree.simulationRole = SimulationRole::InclusiveJet;
  else if (role == "minimum_bias") tree.simulationRole = SimulationRole::MinimumBias;
  else tree.simulationRole = SimulationRole::None;

  tree.outputFile = job.outputFile;

  if (inputMode != "archived_g4_only" && inputMode != "reconstructed_dst") fail("unknown input mode: " + inputMode);
  if (role != "none" && role != "photon_jet" && role != "inclusive_jet" && role != "minimum_bias") fail("unknown simulation role: " + role);

  // ---- source binding -----------------------------------------------------

  tree.source.dataset = job.dataset.empty() ? optionalString(profile, "dataset_token") : job.dataset;
  tree.source.sample = job.sample;
  tree.source.period = job.period;
  tree.source.siDiRole = job.siDiRole;
  tree.source.inputUriSha256 = job.inputUriSha256;
  tree.source.inputFileSha256 = job.inputFileSha256;
  tree.source.sourceManifestSha256 = job.sourceManifestSha256;
  tree.source.run = job.run;
  tree.source.segment = job.segment;
  tree.source.sourceFileOrdinal = job.sourceFileOrdinal;
  tree.source.firstEntry = job.firstEntry;
  tree.source.maxEvents = job.nEvents;

  // ---- nodes --------------------------------------------------------------

  const YAML::Node nodes = root["nodes"];
  tree.nodes.eventHeader = optionalString(nodes, "event_header", "EventHeader");
  tree.nodes.gl1Packet = optionalString(nodes, "gl1_packet", "GL1Packet");
  tree.nodes.globalVertexMap = optionalString(nodes, "global_vertex_map", "GlobalVertexMap");
  tree.nodes.centralityInfo = optionalString(nodes, "centrality_info", "CentralityInfo");
  tree.nodes.minimumBiasInfo = optionalString(nodes, "minimum_bias_info", "MinimumBiasInfo");
  tree.nodes.mbdOut = optionalString(nodes, "mbd_out", "MbdOut");
  tree.nodes.mbdPmts = optionalString(nodes, "mbd_pmts", "MbdPmtContainer");
  tree.nodes.cemcTowers = optionalString(nodes, "cemc_towers", "TOWERINFO_CALIB_CEMC");
  tree.nodes.ihcalTowers = optionalString(nodes, "ihcal_towers", "TOWERINFO_CALIB_HCALIN");
  tree.nodes.ohcalTowers = optionalString(nodes, "ohcal_towers", "TOWERINFO_CALIB_HCALOUT");
  tree.nodes.cemcGeometry = optionalString(nodes, "cemc_geometry", "TOWERGEOM_CEMC");
  tree.nodes.truthInfo = optionalString(nodes, "truth_info", "G4TruthInfo");
  tree.nodes.hepmcEventMap = optionalString(nodes, "hepmc_event_map", "PHHepMCGenEventMap");

  // ---- calorimeter event state ---------------------------------------------

  tree.calorimeter.requireGoodTowers =
      optional<bool>(root["calorimeter"]["event_energy"], "require_good_towers", true);

  // ---- photons ------------------------------------------------------------

  const YAML::Node photons = root["photons"];
  const YAML::Node domain = photons["capture_domain"];
  tree.photon.minEtGeV = required<double>(domain, "min_et_gev", "photons.capture_domain");
  tree.photon.maxEtGeV = required<double>(domain, "max_et_gev", "photons.capture_domain");
  tree.photon.maxAbsEta = required<double>(domain, "max_abs_eta", "photons.capture_domain");
  tree.photon.objectVertexAbsZMaxCm =
      optional<double>(domain, "object_vertex_abs_z_max_cm", tree.photon.objectVertexAbsZMaxCm);
  tree.photon.rawClusterNode = optionalString(photons, "source_cluster_node", "CLUSTERINFO_CEMC");
  tree.photon.photonNode = optionalString(photons, "output_photon_node", "PHOTONCLUSTER_CEMC");
  tree.photon.topoclusterNode = optionalString(photons, "topocluster_node", "TOPOCLUSTER_ALLCALO");

  const std::vector<double> radii = doubleList(photons["isolation"]["radii"]);
  if (!radii.empty())
  {
    tree.photon.isolationRadii = radii;
  }
  tree.photon.truthIsolationSignalRadius =
      optional<double>(photons["isolation"]["truth"], "signal_radius", 0.30);
  tree.photon.truthIsolationMaxEtGeV =
      required<double>(photons["isolation"]["truth"], "max_et_gev", "photons.isolation.truth");
  tree.photon.truthIsolationCoreRadius =
      required<double>(photons["isolation"]["truth"], "self_removal_radius", "photons.isolation.truth");
  tree.photon.timingSampleNs = required<double>(photons["timing"], "sample_ns", "photons.timing");

  const std::vector<std::string> definitions = stringList(photons["shower_shapes"]["definitions"]);
  if (!definitions.empty())
  {
    tree.showerDefinitions = definitions;
  }

  // ---- reconstruction stages -----------------------------------------------

  reco.inputLayout = reco.isArchivedDoubleInteraction ? InputLayout::ArchivedDoubleInteraction
                     : reco.isData                    ? InputLayout::PairedData
                                                      : InputLayout::ReconstructedSimulation;

  reco.truthJetMode = parseTruthJetMode(
      optional<std::string>(profile, "truth_jets_mode", reco.isSimulation ? "auto" : "none"));

  reco.reconstructEventGeometry = !reco.isEmbedded;
  reco.reconstructZdc = reco.isData;
  reco.reconstructCalorimeter = reco.isData || reco.isEmbedded || reco.isArchivedDoubleInteraction;
  reco.runCalorimeterStatusSkimmer = reco.isData;
  reco.setTowerStatus = reco.reconstructCalorimeter &&
                        (!reco.isEmbedded || optional<bool>(profile, "force_embedded_tower_status", false));
  reco.reconstructCentrality = reco.isAuAu && reco.isData &&
                               optional<bool>(profile, "reconstruct_centrality", true);
  reco.runMinimumBiasClassifier = reco.reconstructCentrality &&
                                  optional<bool>(profile, "require_minimum_bias_classifier", true);
  reco.reconstructAuAuBackground = reco.isAuAu;
  reco.reconstructJets = true;
  reco.reconstructPhotons = true;

  const YAML::Node clusterBuilder = root["calorimeter"]["cluster_builder"];
  reco.clusterThresholdGeV =
      required<double>(clusterBuilder, "threshold_energy_gev", "calorimeter.cluster_builder");
  reco.clusterProfile = expandEnvironment(
      required<std::string>(clusterBuilder["profile"], "path", "calorimeter.cluster_builder.profile"));
  reco.simulationHotTowerMapKey =
      optionalString(root["calorimeter"]["tower_status"], "simulation_hot_map_key", "CEMC_hotTowers_status");

  const YAML::Node systemConfig = reco.isPP ? root["pp"] : root["auau"];

  reco.photonMinEtGeV = tree.photon.minEtGeV;
  reco.photonUseVertexCut = optional<bool>(photons["builder"], "use_vertex_cut", false);
  reco.photonVertexCutCm = optional<double>(photons["builder"], "vertex_cut_cm", 60.0);
  reco.showerShapeTowerMinEnergyGeV =
      optional<double>(systemConfig, "shower_shape_tower_min_energy_gev", reco.isPP ? 0.070 : 0.0);

  const YAML::Node isolation = systemConfig["isolation"];
  const std::string isolationMethod = lower(optionalString(isolation, "method", reco.isPP ? "topocluster" : "sub1"));
  reco.useTopoclusterIsolation = reco.isPP && isolationMethod == "topocluster";
  reco.useSubtractedIsolation = reco.isAuAu;
  reco.topoclusterNode = tree.photon.topoclusterNode;

  if (reco.isAuAu)
  {
    const YAML::Node background = root["auau"]["background_subtraction"];
    reco.retowerFractionCut = optional<double>(background["retower_cemc"], "fraction_cut", 0.5);
    reco.firstPassSeedJetD = optional<int>(background["iteration_1"], "seed_jet_D", 3);
    reco.secondPassSeedPtGeV = optional<double>(background["iteration_2"], "seed_pt_min_gev", 7.0);
    reco.flowModulation = required<int>(background, "flow_modulation", "auau.background_subtraction");
    tree.calorimeter.backgroundFlowMode = reco.flowModulation;
  }

  // ---- jets ---------------------------------------------------------------

  tree.jets.truthMatchMaxDeltaR = optional<double>(root["jets"]["matching"], "max_delta_r", 0.30);

  const YAML::Node calibration = root["jets"]["calibration"];
  cond.applyJetEnergyScale = optional<bool>(calibration, "enabled", true);
  cond.jetEnergyScaleCdbKey = optionalString(calibration, "cdb_key", "JES_Calib_Default");
  cond.jetEnergyScalePayloadSha256 =
      lower(required<std::string>(calibration, "expected_payload_sha256", "jets.calibration"));
  cond.jetEnergyScaleUsesEmFraction = optional<bool>(calibration, "use_em_fraction_calibration", false);
  tree.jets.energyScalePayloadSha256 = cond.jetEnergyScalePayloadSha256;
  tree.jets.energyScaleUsesEmFraction = cond.jetEnergyScaleUsesEmFraction;

  const YAML::Node views = systemConfig["jets"]["views"];
  if (!views || !views.IsSequence())
  {
    fail("missing jet view list for the selected collision system");
  }
  for (const YAML::Node& view : views)
  {
    photonjet::JetNodeConfig jet;
    jet.radius = required<double>(view, "radius", "jet view");
    jet.view = parseJetView(required<std::string>(view, "view", "jet view"));
    jet.rawNode = required<std::string>(view, "raw_node", "jet view");
    jet.correctedNode = required<std::string>(view, "corrected_node", "jet view");
    jet.truthNode = optionalString(view, "truth_node");
    jet.inputIdentity = required<std::string>(view, "view", "jet view");
    jet.subtractionIdentity = jet.inputIdentity;
    tree.jets.nodes.push_back(std::move(jet));
  }

  // ---- centrality -----------------------------------------------------------

  tree.centrality.edges = intList(root["auau"]["centrality"]["analysis_bin_edges_percent"]);
  if (tree.centrality.edges.empty())
  {
    tree.centrality.edges = {0, 20, 50, 80};
  }
  tree.centrality.reconstructDataCentrality = reco.reconstructCentrality;

  if (reco.reconstructCentrality)
  {
    const YAML::Node data = root["auau"]["centrality"]["data"];
    tree.centrality.requireOwnRunPayload = optional<bool>(data, "require_own_run_payload", true);
    tree.centrality.calibrationTag = required<std::string>(data, "calibration_set", "auau.centrality.data");
    cond.centralityCalibrationSet = tree.centrality.calibrationTag;
    cond.centralityManifestSha256 =
        lower(required<std::string>(data, "calibration_manifest_sha256", "auau.centrality.data"));
    cond.centralityLocalPayloads = lower(optionalString(data, "source", "cdb")) != "cdb";

    // A run-scoped override table pins the exact payload triplet per run.
    YAML::Node resolved = data;
    if (const YAML::Node runs = root["auau"]["centrality"]["runs"])
    {
      if (const YAML::Node override = runs[std::to_string(job.run)])
      {
        resolved = override;
      }
    }
    cond.centralityPayloadRun = optional<int>(resolved, "payload_run", job.run);
    cond.centralityDivisions = expandEnvironment(optionalString(resolved["divisions"], "path"));
    cond.centralityDivisionsSha256 = lower(optionalString(resolved["divisions"], "sha256"));
    cond.centralityRunScale = expandEnvironment(optionalString(resolved["run_scale"], "path"));
    cond.centralityRunScaleSha256 = lower(optionalString(resolved["run_scale"], "sha256"));
    cond.centralityVertexScale = expandEnvironment(optionalString(resolved["vertex_scale"], "path"));
    cond.centralityVertexScaleSha256 = lower(optionalString(resolved["vertex_scale"], "sha256"));
    cond.centralityFallbackEvidenceSha256 = lower(optionalString(resolved, "fallback_evidence_sha256"));

    tree.centrality.divisionPayloadPath = cond.centralityDivisions;
    tree.centrality.runScalePayloadPath = cond.centralityRunScale;
    tree.centrality.vertexScalePayloadPath = cond.centralityVertexScale;
  }

  // ---- CEMC status and zero suppression, Au+Au data -------------------------

  if (reco.isAuAu && reco.isData)
  {
    const YAML::Node status = root["auau"]["cemc_status"];
    cond.cemcStatusGuard = optional<bool>(status, "validate_paired_input_status", true);
    cond.cemcOriginalCalorimeterInput = expandEnvironment(optionalString(status, "original_calorimeter_input"));
    cond.cemcStatusRecoveryAuthorized = optional<bool>(status, "allow_missing_status_recovery", false);
    cond.cemcStatusRecoveryPayload = expandEnvironment(optionalString(status["recovery_payload"], "path"));
    cond.cemcStatusRecoveryPayloadSha256 = lower(optionalString(status["recovery_payload"], "sha256"));
  }
  cond.cemcZeroSuppressionCrossCalibration =
      reco.isData && optional<bool>(root["calorimeter"]["cemc"], "zero_suppression_cross_calibration", true);

  // ---- archived double interaction -------------------------------------------

  if (reco.isArchivedDoubleInteraction)
  {
    cond.cdbGlobalTag = required<std::string>(profile, "cdb_global_tag", context);
    cond.cdbTimestamp = required<std::uint64_t>(profile, "cdb_timestamp", context);
    const YAML::Node rng = profile["historical_rng_contract"];
    reco.expectedPedestalSequence = required<int>(rng, "expected_pedestal_sequence", context);
    reco.forbidRecoConstsRandomSeed = optional<bool>(rng, "forbid_reco_consts_random_seed", true);
    for (const YAML::Node& seed : rng["seed_sequence"])
    {
      reco.historicalSeeds.push_back(seed.as<unsigned int>());
    }
  }
  else
  {
    cond.cdbTimestamp = job.run > 0 ? static_cast<std::uint64_t>(job.run)
                                    : static_cast<std::uint64_t>(kArchivedDoubleInteractionRun);
    cond.cdbGlobalTag = optionalString(profile, "cdb_global_tag");
  }

  // ---- weights ---------------------------------------------------------------

  if (const YAML::Node reweight = profile["vertex_reweighting"])
  {
    cond.applyVertexReweight = optional<bool>(reweight, "enabled", false);
    cond.vertexReweightFile = expandEnvironment(optionalString(reweight, "file"));
    cond.vertexReweightFileSha256 = lower(optionalString(reweight, "sha256"));
    cond.vertexReweightHistogram = optionalString(reweight, "histogram");
  }
  tree.weights.applyVertexReweight = cond.applyVertexReweight;
  tree.weights.vertexReweightFile = cond.vertexReweightFile;
  tree.weights.vertexReweightFileSha256 = cond.vertexReweightFileSha256;
  tree.weights.vertexReweightHistogram = cond.vertexReweightHistogram;

  // ---- output tables ----------------------------------------------------------

  const YAML::Node output = root["output"];
  tree.output.writePhotonCells = optional<bool>(output, "write_photon_cell_table", true);
  tree.output.writeIsolationConstituents = optional<bool>(output, "write_isolation_constituent_table", true);
  tree.output.writePhotonJetPairs = optional<bool>(output, "write_photon_jet_pair_table", true);

  // ---- retention and provenance ----------------------------------------------

  tree.retainEveryProducerEncounter = optional<bool>(root["capture"], "retain_every_producer_encounter", true);
  tree.requireEmbeddedMinimumBias = optional<bool>(profile, "require_embedded_minimum_bias", false);

  tree.provenance.productionTag = required<std::string>(root["contract"], "name", "contract");
  tree.provenance.softwareRelease = optionalString(root["runtime"], "sphenix_release");
  tree.provenance.producerGitCommit = optionalString(root["runtime"], "producer_git_commit");
  tree.provenance.producerSourceSha256 = optionalString(root["runtime"], "producer_source_sha256");
  tree.provenance.producerLibrarySha256 = optionalString(root["runtime"], "producer_library_sha256");
  tree.provenance.macroSha256 = optionalString(root["runtime"], "macro_sha256");
  tree.provenance.configurationSha256 = plan.configurationSha256;
  tree.provenance.configurationText = plan.configurationText;
  tree.provenance.calibrationManifestSha256 = cond.centralityManifestSha256;

  // ---- physical inputs ----------------------------------------------------------

  plan.inputs = loadInputs(job, reco);

  if (reco.truthJetMode == TruthJetMode::Auto)
  {
    reco.truthJetMode = hasInputStream(plan, "TRUTH_JETS") ? TruthJetMode::DST : TruthJetMode::Build;
  }

  Validate(plan);
  return plan;
}

//
// Validate
//
void Validate(const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  const Conditions& cond = plan.conditions;

  if (plan.inputs.empty()) fail("no input streams were resolved");
  if (plan.job.firstEntry < 0 || plan.job.firstEntry > std::numeric_limits<int>::max() ||
      plan.job.nEvents < -1 || plan.job.nEvents == 0 || plan.job.nEvents > std::numeric_limits<int>::max())
    fail("invalid event range: use -1 for all events, or a positive bounded count");
  // One synchronized source bundle per invocation. Multiple list rows need a
  // source-boundary observer before they can safely share one SourceRecord.
  for (const auto& stream : plan.inputs)
    if (stream.files.size() != 1) fail("one source bundle per invocation is required: " + stream.name);
  if (cond.cdbGlobalTag.empty()) fail("profile must bind an explicit CDB global tag");
  if (!cond.applyJetEnergyScale) fail("canonical reconstructed jets require JES");
  if (cond.jetEnergyScaleCdbKey != "JES_Calib_Default") fail("JES key must match the key consumed by JetCalib");
  if (cond.applyVertexReweight && !validSha256(cond.vertexReweightFileSha256)) fail("vertex weight payload requires SHA256");
  for (const InputStream& stream : plan.inputs)
  {
    if (stream.required && stream.files.empty())
    {
      fail("required input stream is empty: " + stream.name);
    }
  }

  if (reco.isData && plan.tree.simulationRole != SimulationRole::None)
  {
    fail("a data profile cannot carry a simulation role");
  }
  if (reco.isSimulation && plan.tree.simulationRole == SimulationRole::None)
  {
    fail("a simulation profile must name its sample role");
  }

  if (reco.isArchivedDoubleInteraction)
  {
    if (!reco.isPP || !reco.isSimulation) fail("archived double interaction is p+p simulation only");
    // TODO(CODE): port the qualified G4-hit -> waveform/tower sequence and its
    // seed consumers before enabling this lane. CaloTowerCalib alone cannot
    // turn archived G4 hits into towers. Never manufacture an empty success.
    fail("archived DI waveform/tower registration is not yet bound");
    if (reco.truthJetMode != TruthJetMode::DST) fail("archived double interaction requires shipped truth jets");
    if (cond.cdbTimestamp != static_cast<std::uint64_t>(kArchivedDoubleInteractionRun))
    {
      fail("archived double interaction requires the run-28 conditions timestamp");
    }
  }

  if ((reco.isData || reco.isEmbedded) && plan.job.run <= kDataRunThreshold)
  {
    fail("data and embedded production require a real data run number");
  }

  if (reco.reconstructCentrality)
  {
    if (cond.centralityCalibrationSet != kNominalCentralitySet)
    {
      fail("Au+Au centrality is not bound to the nominal frozen calibration set");
    }
    if (cond.centralityManifestSha256 != kNominalCentralityManifestSha256)
    {
      fail("Au+Au centrality manifest SHA256 does not match the frozen nominal set");
    }
    if (cond.centralityLocalPayloads)
    {
      if (cond.centralityPayloadRun <= 0) fail("local centrality requires a positive payload run");
      if (cond.centralityDivisions.empty() || cond.centralityRunScale.empty() || cond.centralityVertexScale.empty())
      {
        fail("local centrality requires divisions, run-scale and vertex-scale payloads");
      }
      if (!validSha256(cond.centralityDivisionsSha256) || !validSha256(cond.centralityRunScaleSha256) ||
          !validSha256(cond.centralityVertexScaleSha256))
      {
        fail("local centrality requires an exact SHA256 for each of the three payloads");
      }
      if (cond.centralityPayloadRun != plan.job.run)
      {
        if (cond.centralityPayloadRun != kQualifiedCentralityFallbackRun)
        {
          fail("only the qualified run-68144 centrality fallback is permitted");
        }
        if (!validSha256(cond.centralityFallbackEvidenceSha256))
        {
          fail("the centrality fallback requires pinned qualification evidence");
        }
      }
      else if (!cond.centralityFallbackEvidenceSha256.empty())
      {
        fail("own-run centrality cannot carry a fallback approval");
      }
    }
  }

  if (cond.cemcStatusGuard)
  {
    if (cond.cemcOriginalCalorimeterInput.empty())
    {
      fail("the CEMC status guard requires the exact original calorimeter input file");
    }
    if (cond.cemcStatusRecoveryAuthorized &&
        (cond.cemcStatusRecoveryPayload.empty() || !validSha256(cond.cemcStatusRecoveryPayloadSha256)))
    {
      fail("CEMC status recovery requires a pinned run-specific payload and its SHA256");
    }
  }

  if (cond.applyJetEnergyScale)
  {
    if (cond.jetEnergyScaleUsesEmFraction)
    {
      fail("the nominal jet energy scale is the legacy method without the electromagnetic-fraction calibration");
    }
    if (cond.jetEnergyScalePayloadSha256 != kNominalJetEnergyScaleSha256)
    {
      fail("jet energy-scale payload SHA256 differs from the nominal common p+p v6 payload");
    }
  }

  if (cond.applyVertexReweight)
  {
    if (!reco.isSimulation) fail("vertex reweighting applies to simulation only");
    if (cond.vertexReweightFile.empty() || cond.vertexReweightHistogram.empty())
    {
      fail("vertex reweighting requires a file and a histogram name");
    }
  }

  if (reco.isAuAu && (reco.flowModulation < 0 || reco.flowModulation > 3))
  {
    fail("Au+Au flow modulation must be resolved to the accepted production value");
  }
}

//
// ConfigureConditions
//
void ConfigureConditions(Plan& plan)
{
  Conditions& cond = plan.conditions;
  const Reconstruction& reco = plan.reconstruction;

  recoConsts* rc = recoConsts::instance();
  const int run = reco.isArchivedDoubleInteraction ? kArchivedDoubleInteractionRun : plan.job.run;
  rc->set_IntFlag("RUNNUMBER", run);
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
}

//
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
  registerPhotons(server, plan);

  // 9. The producer, last, so it reads a complete event.
  server->registerSubsystem(producer);

  (void) reco;
}

//
// RegisterInputs
//
void RegisterInputs(Fun4AllServer* server, const Plan& plan)
{
  if (!server) fail("null Fun4AllServer");

  // Reconstructing the calorimeter needs the tower geometry run node.
  if (plan.reconstruction.reconstructCalorimeter)
  {
    auto* geometry = new Fun4AllRunNodeInputManager("DST_GEO");
    geometry->AddFile(CDBInterface::instance()->getUrl("calo_geo"));
    server->registerInputManager(geometry);
  }

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
    Fun4AllInputManager* manager = makeInputManager(stream);
    manager->Verbosity(plan.job.verbosity);
    for (const std::string& file : stream.files)
    {
      manager->AddFile(file);
    }
    server->registerInputManager(manager);
  }
}

//
// Print
//
void Print(const Plan& plan)
{
  const Reconstruction& reco = plan.reconstruction;
  const Conditions& cond = plan.conditions;

  const char* truthMode = "none";
  switch (reco.truthJetMode)
  {
    case TruthJetMode::None: truthMode = "none"; break;
    case TruthJetMode::Auto: truthMode = "auto"; break;
    case TruthJetMode::DST: truthMode = "dst"; break;
    case TruthJetMode::Build: truthMode = "build"; break;
  }

  std::cout << "\n"
            << "============================================================\n"
            << " PhotonJetTree production plan\n"
            << "============================================================\n"
            << " profile              : " << plan.job.profile << "\n"
            << " system               : " << (reco.isPP ? "pp" : "AuAu") << "\n"
            << " data kind            : " << (reco.isData ? "DATA" : "SIM") << "\n"
            << " embedded             : " << (reco.isEmbedded ? "yes" : "no") << "\n"
            << " archived DI          : " << (reco.isArchivedDoubleInteraction ? "yes" : "no") << "\n"
            << " run / segment        : " << plan.job.run << " / " << plan.job.segment << "\n"
            << " first entry / events : " << plan.job.firstEntry << " / " << plan.job.nEvents << "\n"
            << " input list           : " << plan.job.inputList << "\n"
            << " output               : " << plan.job.outputFile << "\n"
            << " config sha256        : " << plan.configurationSha256 << "\n"
            << " truth jets           : " << truthMode << "\n"
            << " jet collections      : " << plan.tree.jets.nodes.size() << "\n"
            << " CDB timestamp / tag  : " << cond.cdbTimestamp << " / "
            << (cond.cdbGlobalTag.empty() ? "<preconfigured>" : cond.cdbGlobalTag) << "\n"
            << " calorimeter          : " << (reco.reconstructCalorimeter ? "reconstructed here" : "from input") << "\n"
            << " JES                  : " << (cond.applyJetEnergyScale ? "legacy, once, " + cond.jetEnergyScalePayload : "disabled") << "\n";

  if (reco.reconstructCentrality)
  {
    std::cout << " centrality           : " << cond.centralityCalibrationSet
              << (cond.centralityLocalPayloads ? " local payload run " + std::to_string(cond.centralityPayloadRun)
                                               : " from CDB")
              << "\n";
  }
  if (cond.cemcStatusGuard)
  {
    std::cout << " CEMC status guard    : " << cond.cemcOriginalCalorimeterInput
              << (cond.cemcStatusRecoveryAuthorized ? " (recovery authorised)" : "") << "\n";
  }
  if (cond.applyVertexReweight)
  {
    std::cout << " vertex reweight      : " << cond.vertexReweightFile << ":" << cond.vertexReweightHistogram << "\n";
  }

  std::cout << " input streams        : " << plan.inputs.size() << "\n";
  for (const InputStream& stream : plan.inputs)
  {
    std::cout << "   " << std::left << std::setw(12) << stream.name
              << " files=" << stream.files.size()
              << " sync=" << (stream.synchronized ? "yes" : "no") << "\n";
  }
  std::cout << "============================================================\n" << std::endl;
}

}  // namespace production
}  // namespace photonjet
