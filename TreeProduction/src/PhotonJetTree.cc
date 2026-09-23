//
// PhotonJetTree.cc
//
// Lifecycle and event flow of the isolated-photon + jet tree producer.
//
// This file owns the order in which things happen and the way failures are
// reported. Every physics step it calls lives in one of the internal files:
//
//   internal/Event.cc     identity, trigger, vertex, minimum bias, centrality,
//                         calorimeter state, weights, event accounting
//   internal/Photons.cc   photon candidates, timing, shower shapes, cells,
//                         isolation
//   internal/Jets.cc      reconstructed jets and photon-jet pairs
//   internal/Truth.cc     simulation truth and reconstruction-to-truth links
//   internal/Output.cc    ROOT tables, serialization, provenance
//
// Failure policy
//
//   A configuration or setup problem stops the run before any event is read.
//   A contract violation during an event, for example a required node that is
//   absent in a mode that needs it, also stops the run: silently skipping the
//   event would corrupt the exposure accounting the trees are built to
//   preserve. Ordinary physics absence (no photon candidate, no truth jet) is
//   not a failure and produces a normal event row with empty object tables.
//
#include "PhotonJetTree.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <exception>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <TH1.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{

// Report a failure at a named lifecycle stage and return the Fun4All code that
// stops the run. The message names the module so a log with several modules
// stays readable.
int abortRun(const std::string& module,
             const char* stage,
             const std::exception& error)
{
  std::cerr << "\n"
            << module << ": " << stage << " failed\n"
            << "  " << error.what() << "\n"
            << "  the run is stopped; no partial output is certified\n"
            << std::endl;

  return Fun4AllReturnCodes::ABORTRUN;
}

}  // namespace

//
// Construction
//
PhotonJetTree::PhotonJetTree(Config config, const std::string& name)
  : SubsysReco(name)
  , m_config(std::move(config))
{
}

PhotonJetTree::~PhotonJetTree()
{
  // End() is the supported way to finish. If it never ran, the file is closed
  // without completion metadata, which downstream readers treat as an
  // incomplete part rather than as a short one.
  if (m_outputFile && !m_finalized)
  {
    std::cerr << Name()
              << ": destroyed before End(); output is not marked complete"
              << std::endl;

    try { closeOutput(); }
    catch (const std::exception& error) { std::cerr << error.what() << std::endl; }
  }
  delete m_vertexReweight;
}

//
// Production mode
//
bool PhotonJetTree::isPP() const
{
  return m_config.system == CollisionSystem::PP;
}

bool PhotonJetTree::isAuAu() const
{
  return m_config.system == CollisionSystem::AuAu;
}

bool PhotonJetTree::isData() const
{
  return m_config.dataKind == DataKind::Data;
}

bool PhotonJetTree::isSimulation() const
{
  return m_config.dataKind == DataKind::Simulation;
}

// Au+Au simulation in this analysis is always a hard-scatter event embedded in
// a real minimum-bias event; there is no stand-alone Au+Au generator sample.
bool PhotonJetTree::isEmbeddedSimulation() const
{
  return isAuAu() && isSimulation();
}

bool PhotonJetTree::isArchivedDoubleInteraction() const
{
  return m_config.inputMode == InputMode::ArchivedG4Only;
}

bool PhotonJetTree::isPhotonJetSimulation() const
{
  return isSimulation() &&
         m_config.simulationRole == SimulationRole::PhotonJet;
}

bool PhotonJetTree::isInclusiveJetSimulation() const
{
  return isSimulation() &&
         m_config.simulationRole == SimulationRole::InclusiveJet;
}

//
// Init
//
//   validate configuration
//     -> initialize source record and identity
//     -> open output, book tables, write file provenance
//
int PhotonJetTree::Init(PHCompositeNode* /*topNode*/)
{
  try
  {
    validateConfiguration();
    initializeSourceRecord();
    loadWeightInputs();
    initializeOutput();

    m_outputReady = true;
  }
  catch (const std::exception& error)
  {
    markAborted();
    return abortRun(Name(), "Init", error);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//
// InitRun
//
//   capture run-dependent state (trigger configuration table)
//
int PhotonJetTree::InitRun(PHCompositeNode* topNode)
{
  if (!m_outputReady)
  {
    std::cerr << Name() << ": InitRun called before Init succeeded" << std::endl;
    markAborted();
    return Fun4AllReturnCodes::ABORTRUN;
  }

  try
  {
    initializeRun(topNode);

    m_runInitialized = true;
  }
  catch (const std::exception& error)
  {
    markAborted();
    return abortRun(Name(), "InitRun", error);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//
// process_event
//
// The physical order of the analysis, written out in full.
//
int PhotonJetTree::process_event(PHCompositeNode* topNode)
{
  if (!m_outputReady || !m_runInitialized)
  {
    std::cerr << Name() << ": event received before setup completed" << std::endl;
    markAborted();
    return Fun4AllReturnCodes::ABORTRUN;
  }

  try
  {
    // Every event this module is handed is an encounter, whatever it contains.
    // The ordinal is assigned before any capture so that a failure later in
    // the event still leaves the accounting consistent. The pending upstream
    // observation, if the exposure observer is registered, is consumed here.
    resetEventState();

    m_event.producerEventOrdinal = m_source.encounteredEvents;
    ++m_source.encounteredEvents;

    m_upstreamPending = false;

    //
    // Event context
    //
    captureEventIdentity(topNode);
    captureTriggerInformation(topNode);
    captureVertexInformation(topNode);
    captureMbdInformation(topNode);
    captureCentrality(topNode);
    captureCalorimeterInformation(topNode);

    //
    // Reconstructed objects
    //
    capturePhotons(topNode);
    captureJets(topNode);
    buildPhotonJetPairs();

    //
    // Simulation truth and its relation to reconstruction
    //
    // Applicability is a property of the production mode. Data leaves the
    // truth tables empty by construction, not because a node was missing.
    //
    if (isSimulation())
    {
      captureSimulationTruth(topNode);
      buildPhotonTruthAssociations(topNode);
      buildJetTruthAssociations();
    }

    //
    // Producer weight
    //
    captureEventWeights(topNode);

    //
    // Output
    //
    // Every encounter is an exposure record, including zero-candidate events.
    fillEventOutput();
    fillObjectOutput();
    ++m_source.retainedEvents;
  }
  catch (const std::exception& error)
  {
    markAborted();
    return abortRun(Name(), "process_event", error);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//
// ResetEvent
//
// Only per-event state is cleared. Source accounting, scaler snapshots, the
// trigger configuration table and the output handles persist.
//
int PhotonJetTree::ResetEvent(PHCompositeNode* /*topNode*/)
{
  resetEventState();

  return Fun4AllReturnCodes::EVENT_OK;
}

//
// End
//
//   flush remaining source-level rows
//     -> write the source record and completion metadata
//     -> close the file
//
int PhotonJetTree::End(PHCompositeNode* /*topNode*/)
{
  if (m_finalized)
  {
    return m_aborted ? Fun4AllReturnCodes::ABORTRUN : Fun4AllReturnCodes::EVENT_OK;
  }

  if (!m_outputReady)
  {
    std::cerr << Name() << ": End called without a successful Init" << std::endl;
    markAborted();
    return Fun4AllReturnCodes::ABORTRUN;
  }

  try
  {
    // An observation still pending after the last event belongs to an event
    // that reconstruction stopped; account for it before the source row is
    // written.
    finishUpstreamAccounting();

    writeCompletionMetadata();
    closeOutput();

    m_finalized = true;
  }
  catch (const std::exception& error)
  {
    markAborted();
    return abortRun(Name(), "End", error);
  }

  std::cout << Name()
            << ": encountered " << m_source.encounteredEvents
            << ", retained " << m_source.retainedEvents
            << ", upstream rejected " << m_source.upstreamRejectedEvents
            << ", output " << m_config.outputFile
            << std::endl;

  return m_aborted ? Fun4AllReturnCodes::ABORTRUN : Fun4AllReturnCodes::EVENT_OK;
}

//
// Upstream exposure accounting
//
void PhotonJetTree::recordUpstreamRejectedEvent(
    const int run,
    const std::int64_t physicalEventSequence,
    const std::int64_t sourceEntry)
{
  UpstreamRejectedEventRecord rejected;

  rejected.sourceId = m_source.sourceId;
  rejected.run = run;
  rejected.physicalEventSequence = physicalEventSequence;
  rejected.sourceEntry = sourceEntry;

  m_upstreamRejectedEvents.push_back(rejected);

  ++m_source.upstreamRejectedEvents;
}

//
// Configuration validation
//
// Everything checked here can be checked before a file is opened. The first
// inconsistency found is reported; there is no value in listing several when
// the job cannot run.
//
void PhotonJetTree::validateConfiguration() const
{
  const Config& c = m_config;

  if (c.outputFile.empty())
  {
    throw std::runtime_error("output file name is empty");
  }

  // Mode consistency.

  if (isData() && c.simulationRole != SimulationRole::None)
  {
    throw std::runtime_error("a data production cannot carry a simulation role");
  }

  if (isSimulation() && c.simulationRole == SimulationRole::None)
  {
    throw std::runtime_error("a simulation production must name its sample role");
  }

  if (isArchivedDoubleInteraction() && !(isPP() && isSimulation()))
  {
    throw std::runtime_error(
        "archived Geant4-only input is only defined for p+p simulation");
  }

  if (c.requireEmbeddedMinimumBias && !isEmbeddedSimulation())
  {
    throw std::runtime_error(
        "the embedded minimum-bias requirement applies only to embedded simulation");
  }

  // Source identity needs a pinned digest. Without it the compact keys have no
  // provenance and the output could not be joined with anything.

  auto digestValid = [](const std::string& value) {
    return value.size() == 64 && std::all_of(value.begin(), value.end(),
        [](unsigned char c) { return std::isxdigit(c); });
  };
  for (const auto* digest : {&c.provenance.producerSourceSha256,
                             &c.provenance.producerLibrarySha256,
                             &c.provenance.macroSha256,
                             &c.provenance.configurationSha256})
    if (!digestValid(*digest)) throw std::runtime_error("source/library/macro/configuration provenance is not pinned");
  const bool haveFileDigest = digestValid(c.source.inputFileSha256);
  const bool haveManifestDigest = digestValid(c.source.sourceManifestSha256);
  if ((!c.source.inputFileSha256.empty() && !haveFileDigest) ||
      (!c.source.sourceManifestSha256.empty() && !haveManifestDigest))
    throw std::runtime_error("supplied source digest is malformed");
  if (!haveFileDigest && haveManifestDigest && c.source.sourceFileOrdinal < 0)
    throw std::runtime_error("manifest-based source identity requires a source ordinal");
  if (!c.retainEveryProducerEncounter)
    throw std::runtime_error("canonical trees must retain every producer encounter");
  for (double value : {c.photon.minEtGeV, c.photon.maxEtGeV, c.photon.maxAbsEta,
                      c.photon.objectVertexAbsZMaxCm, c.jets.truthMatchMaxDeltaR,
                      c.photon.truthIsolationSignalRadius, c.photon.truthIsolationMaxEtGeV,
                      c.photon.truthIsolationCoreRadius})
    if (!std::isfinite(value) || value < 0)
      throw std::runtime_error("storage and matching bounds must be finite and nonnegative");
  if (std::abs(c.photon.truthIsolationSignalRadius - 0.30) > 1e-6 &&
      std::abs(c.photon.truthIsolationSignalRadius - 0.40) > 1e-6)
    throw std::runtime_error("truth isolation is stored at R=0.3 and R=0.4 only");
  std::set<std::pair<JetView, std::uint64_t>> collections;
  for (const auto& jet : c.jets.nodes)
    if (!std::isfinite(jet.radius) || jet.radius <= 0 ||
        !collections.insert({jet.view, photonjet::radiusCode(jet.radius)}).second)
      throw std::runtime_error("invalid or duplicate jet view/radius");

  if (!haveFileDigest && !haveManifestDigest)
  {
    throw std::runtime_error(
        "source identity requires a pinned input-file or source-manifest SHA256");
  }

  if (c.source.firstEntry < 0)
  {
    throw std::runtime_error("first source entry cannot be negative");
  }

  // Photon storage domain.

  if (!(c.photon.minEtGeV < c.photon.maxEtGeV))
  {
    throw std::runtime_error("photon capture domain has minEt >= maxEt");
  }

  if (!(c.photon.maxAbsEta > 0.0))
  {
    throw std::runtime_error("photon capture domain has non-positive maxAbsEta");
  }

  if (c.photon.isolationRadii.empty())
  {
    throw std::runtime_error("at least one isolation radius must be configured");
  }

  std::set<std::uint64_t> isolationCodes;
  for (const double radius : c.photon.isolationRadii)
  {
    if (!std::isfinite(radius) || !isolationCodes.insert(photonjet::radiusCode(radius)).second || (std::abs(radius - 0.3) > 1e-6 && std::abs(radius - 0.4) > 1e-6))
    {
      throw std::runtime_error("only R=0.3 and R=0.4 isolation cones are supported");
    }
  }

  if (!std::isfinite(c.photon.timingSampleNs) || !(c.photon.timingSampleNs > 0.0))
  {
    throw std::runtime_error("timing sample length must be positive");
  }

  if (!(c.photon.truthIsolationSignalRadius > 0.0))
  {
    throw std::runtime_error("truth isolation signal radius must be positive");
  }

  // Jet collections.

  if (c.jets.nodes.empty())
  {
    throw std::runtime_error("no reconstructed jet collection is configured");
  }

  for (const JetNodeConfig& jet : c.jets.nodes)
  {
    if (jet.rawNode.empty() || jet.correctedNode.empty())
    {
      throw std::runtime_error("a jet collection has an empty raw or corrected node");
    }

    if (jet.rawNode == jet.correctedNode)
    {
      throw std::runtime_error(
          "a jet collection names the same node before and after the energy scale");
    }

    if (!(jet.radius > 0.0))
    {
      throw std::runtime_error("a jet collection has a non-positive radius");
    }

    const bool ppView = jet.view == JetView::PP;

    if (ppView != isPP())
    {
      throw std::runtime_error(
          "a jet collection view does not belong to the configured collision system");
    }
  }

  if (!(c.jets.truthMatchMaxDeltaR > 0.0))
  {
    throw std::runtime_error("truth jet matching radius must be positive");
  }

  // Centrality.

  if (c.centrality.reconstructDataCentrality && !(isAuAu() && isData()))
  {
    throw std::runtime_error(
        "centrality reconstruction during this pass applies only to Au+Au data");
  }

  for (double edge : c.centrality.edges)
    if (!std::isfinite(edge)) throw std::runtime_error("centrality edges must be finite");
  for (std::size_t i = 1; i < c.centrality.edges.size(); ++i)
  {
    if (c.centrality.edges[i] <= c.centrality.edges[i - 1])
    {
      throw std::runtime_error("centrality bin edges must increase strictly");
    }
  }
}

//
// Source record
//
void PhotonJetTree::initializeSourceRecord()
{
  const SourceConfig& s = m_config.source;

  m_source = SourceRecord{};

  m_source.dataset = s.dataset;
  m_source.sample = s.sample;
  m_source.period = s.period;
  m_source.siDiRole = s.siDiRole;

  m_source.inputUriSha256 = s.inputUriSha256;
  m_source.inputFileSha256 = s.inputFileSha256;
  m_source.sourceManifestSha256 = s.sourceManifestSha256;

  m_source.run = s.run;
  m_source.segment = s.segment;
  m_source.sourceFileOrdinal = s.sourceFileOrdinal;

  m_source.firstEntry = s.firstEntry;
  m_source.lastEntry = -1;

  // The file digest names the source when there is one. A manifest digest is
  // the fallback for sources that are not one physical file.
  const std::string& digest =
      s.inputFileSha256.size() >= 32 ? s.inputFileSha256
                                     : s.sourceManifestSha256;

  m_source.sourceId = photonjet::makeSourceIdentity(digest);
  if (s.inputFileSha256.empty())
  {
    m_source.sourceId.hi = photonjet::mixIdentityWord(m_source.sourceId.hi ^ static_cast<std::uint64_t>(s.sourceFileOrdinal));
    m_source.sourceId.lo = photonjet::mixIdentityWord(m_source.sourceId.lo ^ static_cast<std::uint64_t>(s.sourceFileOrdinal));
  }

  if (!m_source.sourceId.valid())
  {
    throw std::runtime_error(
        "could not build a source identity from the configured SHA256");
  }
}

//
// Run-dependent state
//
void PhotonJetTree::initializeRun(PHCompositeNode* topNode)
{
  // The trigger configuration is a property of the run, not of an event, and
  // simulation has none.
  if (isData())
  {
    captureTriggerRunInformation(topNode);
  }
}
