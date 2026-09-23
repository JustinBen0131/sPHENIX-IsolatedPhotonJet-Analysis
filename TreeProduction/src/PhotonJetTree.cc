/**
 * @file PhotonJetTree.cc
 * @brief Lifecycle and event-flow spine of the isolated-photon + jet producer.
 *
 * PhotonJetTree is one Fun4All subsystem serving pp and AuAu, DATA and simulation. This file controls the lifecycle and the fixed order in which
 * one event is converted into canonical tree records; the physics algorithms themselves live in the responsibility-specific implementation files.
 *
 * Event flow
 * ----------
 *   reset
 *     -> source / event identity
 *     -> trigger and scaler state
 *     -> reconstructed vertex and MBD state
 *     -> centrality
 *     -> calorimeter event state
 *     -> reconstructed photons and isolation
 *     -> reconstructed jets and photon-jet pairs
 *     -> simulation truth and reco-truth associations
 *     -> producer weights
 *     -> ROOT output
 *
 * Implementation ownership
 * ------------------------
 *   internal/Event.cc
 *       Event identity, triggers/scalers, vertex, MBD, centrality,
 *       calorimeter state, weights, and exposure accounting.
 *
 *   internal/Photons.cc
 *       Photon candidates, timing, shower-shape views, cells, and isolation.
 *
 *   internal/Jets.cc
 *       Reconstructed jet collections and photon-jet relationships.
 *
 *   internal/Truth.cc
 *       Simulation truth and reconstruction-to-truth associations.
 *
 *   internal/Output.cc
 *       ROOT booking, serialization, provenance, and completion metadata.
 *
 * Failure policy
 * --------------
 * Setup or production-contract failures abort the run rather than silently changing the event population. Ordinary physics absence is not an error:
 * an event with no photon, jet, or truth match is still a valid event record with the corresponding object collections empty.
 *
 * This file deliberately does not own reconstruction algorithms, photon-ID models, working points, analysis selections, histogramming, or unfolding.
 */


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
/**
 * Convert a lifecycle exception into a consistent Fun4All run abort.
 *
 * Keeping failure reporting here gives Init(), InitRun(), process_event(), and
 * End() one common policy: a broken production contract stops the run and the
 * output is never presented as a completed product.
 */
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


// ============================================================================
// Construction and ownership
// ============================================================================

PhotonJetTree::PhotonJetTree(Config config, const std::string& name)
  : SubsysReco(name)
  , m_config(std::move(config))
{
}

PhotonJetTree::~PhotonJetTree()
{
  /*
   * End() is the only path that certifies a completed product.
   *
   * If destruction happens earlier, close the physical file but deliberately
   * omit completion metadata. Downstream consumers can therefore distinguish
   * an interrupted job from a legitimately short production.
   */
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


// ============================================================================
// Production-mode queries
//
// These helpers keep mode-dependent logic readable throughout the split
// implementation files. They describe the configured production lane; they do
// not infer mode from the accidental presence or absence of DST nodes.
// ============================================================================

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

/*
 * Au+Au simulation in this analysis is embedded simulation: a simulated hard
 * scattering is overlaid with a minimum-bias Au+Au event. There is no
 * stand-alone Au+Au generator production lane represented by this mode.
 */
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


// ============================================================================
// Fun4All lifecycle
// ============================================================================

/**
 * Initialize immutable job state and the output product.
 *
 * Order:
 *   configuration validation
 *     -> source identity / provenance
 *     -> external weight inputs
 *     -> ROOT file and table booking
 *
 * No event is read until this sequence completes successfully.
 */
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


/**
 * Initialize state whose meaning is tied to the current run.
 *
 * At present this primarily captures the DATA trigger configuration. The
 * method remains a distinct lifecycle boundary so run-dependent state never
 * leaks into ordinary per-event capture.
 */
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


/**
 * Capture one producer encounter.
 *
 * The sequence below is intentionally explicit: this function is the readable
 * spine of the DST -> canonical-tree transformation.
 *
 *   accounting / identity
 *     -> event context
 *     -> reconstructed photons
 *     -> reconstructed jets and recoil relations
 *     -> simulation truth and associations
 *     -> producer weights
 *     -> serialization
 *
 * DATA and simulation share this path. Truth work is skipped by declared mode,
 * not by treating missing truth nodes as an implicit DATA signal.
 */
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
    /*
     * Start every producer encounter with a fresh event buffer and assign its
     * ordinal before any physics capture. Exposure accounting therefore does
     * not depend on whether this event ultimately contains photons or jets.
     *
     * Clearing m_upstreamPending also consumes the observer handoff for the
     * event that successfully reached this producer.
     */
    resetEventState();
    m_event.producerEventOrdinal = m_source.encounteredEvents;
    ++m_source.encounteredEvents;
    m_upstreamPending = false;

    // ------------------------------------------------------------------------
    // Event context
    //
    // Establish the event identity first; every object captured below is keyed
    // to this context.
    // ------------------------------------------------------------------------
    captureEventIdentity(topNode);
    captureTriggerInformation(topNode);
    captureVertexInformation(topNode);
    captureMbdInformation(topNode);
    captureCentrality(topNode);
    captureCalorimeterInformation(topNode);

    // ------------------------------------------------------------------------
    // Reconstructed objects
    //
    // Photons and jets are captured independently; pair records are formed only
    // after both reconstructed collections exist.
    // ------------------------------------------------------------------------
    capturePhotons(topNode);
    captureJets(topNode);
    buildPhotonJetPairs();

    // ------------------------------------------------------------------------
    // Simulation truth and reconstruction-to-truth relations
    //
    // Applicability is declared by the production mode. DATA therefore has
    // intentionally empty truth collections rather than "missing" truth.
    // ------------------------------------------------------------------------
    if (isSimulation())
    {
      captureSimulationTruth(topNode);
      buildPhotonTruthAssociations(topNode);
      buildJetTruthAssociations();
    }

    // ------------------------------------------------------------------------
    // Event weight
    //
    // Capture the producer-owned event weight after event/truth context is
    // available. Downstream sample stitching and final normalization remain
    // outside TreeProduction.
    // ------------------------------------------------------------------------
    captureEventWeights(topNode);

    // ------------------------------------------------------------------------
    // Persist the complete encounter
    //
    // The event row is written even when every object collection is empty.
    // This preserves the producer-level exposure denominator.
    // ------------------------------------------------------------------------
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


/**
 * Clear per-event buffers after an event.
 *
 * Source accounting, run-level trigger/scaler state, output handles, and other
 * job-lifetime state intentionally survive this reset.
 */
int PhotonJetTree::ResetEvent(PHCompositeNode* /*topNode*/)
{
  resetEventState();
  return Fun4AllReturnCodes::EVENT_OK;
}


/**
 * Finalize and certify the output product.
 *
 * Remaining upstream accounting is resolved before completion metadata is
 * written. Completion metadata is therefore the terminal declaration that the
 * source accounting and ROOT serialization finished successfully.
 */
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
    /*
     * A still-pending upstream observation represents an input event that did
     * not reach process_event(). Resolve it before the terminal source counts
     * are serialized.
     */
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


// ============================================================================
// Upstream exposure accounting
// ============================================================================

/**
 * Record an input event rejected before PhotonJetTree::process_event().
 *
 * These rows preserve the distinction between the original source population
 * and the subset that reached this producer. They are exposure accounting, not
 * reconstructed-event records.
 */
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


// ============================================================================
// Configuration validation
// ============================================================================

/**
 * Validate the complete producer contract before event processing begins.
 *
 * This function checks only configuration invariants that can be established
 * without reading an event: mode consistency, provenance, storage domains,
 * supported isolation radii, jet collections, and centrality configuration.
 *
 * Validation is deliberately fail-fast. A contradictory configuration must
 * never become a partially interpretable ROOT product.
 */
void PhotonJetTree::validateConfiguration() const
{
  const Config& c = m_config;

  // --------------------------------------------------------------------------
  // Output and production-mode consistency
  // --------------------------------------------------------------------------
  if (c.outputFile.empty())
  {
    throw std::runtime_error("output file name is empty");
  }

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

  // --------------------------------------------------------------------------
  // Provenance and source identity
  //
  // Compact relational identities are only meaningful when their originating
  // source/configuration/software bytes are themselves pinned.
  // --------------------------------------------------------------------------
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

  /*
   * One event row per producer encounter is a core schema invariant. Turning
   * this off would change the event/exposure denominator downstream.
   */
  if (!c.retainEveryProducerEncounter)
    throw std::runtime_error("canonical trees must retain every producer encounter");

  // --------------------------------------------------------------------------
  // Common numeric storage and matching bounds
  // --------------------------------------------------------------------------
  for (double value : {c.photon.minEtGeV, c.photon.maxEtGeV, c.photon.maxAbsEta,
                       c.photon.objectVertexAbsZMaxCm, c.jets.truthMatchMaxDeltaR,
                       c.photon.truthIsolationSignalRadius, c.photon.truthIsolationMaxEtGeV,
                       c.photon.truthIsolationCoreRadius})
    if (!std::isfinite(value) || value < 0)
      throw std::runtime_error("storage and matching bounds must be finite and nonnegative");

  if (std::abs(c.photon.truthIsolationSignalRadius - 0.30) > 1e-6 &&
      std::abs(c.photon.truthIsolationSignalRadius - 0.40) > 1e-6)
    throw std::runtime_error("truth isolation is stored at R=0.3 and R=0.4 only");

  /*
   * A reconstructed jet collection is uniquely identified by its physical view
   * and radius. Duplicate definitions would make downstream joins ambiguous.
   */
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

  // --------------------------------------------------------------------------
  // Photon capture domain
  // --------------------------------------------------------------------------
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

  /*
   * The canonical tree stores the two explicitly supported isolation cones.
   * Working points remain downstream; this only fixes which reconstructed
   * observables are available for later analysis.
   */
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

  // --------------------------------------------------------------------------
  // Reconstructed jet collections
  // --------------------------------------------------------------------------
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

    /*
     * Raw and calibrated collections must remain distinct so the applied JES
     * can be represented explicitly rather than inferred from one mutable node.
     */
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

  // --------------------------------------------------------------------------
  // Centrality
  //
  // Active centrality reconstruction belongs only to Au+Au DATA in this
  // production contract. Other modes either consume their existing state or
  // have no centrality concept.
  // --------------------------------------------------------------------------
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


// ============================================================================
// Source identity and source-level provenance
// ============================================================================

/**
 * Initialize the one source record represented by this producer invocation.
 *
 * A physical file is identified directly by its SHA256. For manifest-addressed
 * sources, the manifest hash plus source ordinal distinguishes the member.
 * The resulting 128-bit source identity is a compact join key; the full hashes
 * remain stored independently as provenance.
 */
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

  /*
   * Prefer the physical file digest whenever one exists. A source represented
   * only through a manifest is anchored to that manifest and disambiguated by
   * its member ordinal below.
   */
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


// ============================================================================
// Run-dependent state
// ============================================================================

/**
 * Capture state that belongs to the run rather than to an individual event.
 */
void PhotonJetTree::initializeRun(PHCompositeNode* topNode)
{
  // Trigger configuration is meaningful for DATA runs and intentionally absent
  // from simulation rather than represented as a failed lookup.
  if (isData())
  {
    captureTriggerRunInformation(topNode);
  }
}
