#ifndef PHOTONJETTREE_H
#define PHOTONJETTREE_H

/**
 * @file PhotonJetTree.h
 * Public producer lifecycle and retained record state.
 * Event, object, truth and output operations are implemented in internal/.
 * Types.h owns record types; the producer never selects candidates by BDT.
 * Empty object collections remain valid retained events.
 */

#include "internal/Types.h"

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class PHCompositeNode;

class TFile;
class TH1;
class TH2;
class TTree;

class RawCluster;
class Jet;
class TowerInfoContainer;
namespace photonjet { class PhotonClusterBuilder; class DirectEventWriter; }
namespace photonjet::production { class PinnedModel; }

class PHG4Particle;
class CaloRawClusterEval;

class PhotonJetTree : public SubsysReco
{
 public:
  // ==========================================================================
  // Canonical producer vocabulary
  //
  // internal/Types.h is the sole owner of these definitions. The aliases below
  // expose that vocabulary through PhotonJetTree without creating a second type
  // system.
  // ==========================================================================

  using CollisionSystem = photonjet::CollisionSystem;
  using DataKind = photonjet::DataKind;
  using SimulationRole = photonjet::SimulationRole;
  using InputMode = photonjet::InputMode;
  using JetView = photonjet::JetView;

  using IsolationMethod = photonjet::IsolationMethod;
  using TriState = photonjet::TriState;
  using CaptureState = photonjet::CaptureState;
  using AssociationState = photonjet::AssociationState;
  using DominantTruthState = photonjet::DominantTruthState;
  using DominantTruthEvaluator = photonjet::DominantTruthEvaluator;
  using TruthPhotonClass = photonjet::TruthPhotonClass;
  using LinkKind = photonjet::LinkKind;

  using Identity = photonjet::Identity;

  using Config = photonjet::Config;
  using WeightConfig = photonjet::WeightConfig;
  using SourceConfig = photonjet::SourceConfig;
  using NodeConfig = photonjet::NodeConfig;
  using CalorimeterConfig = photonjet::CalorimeterConfig;
  using PhotonConfig = photonjet::PhotonConfig;
  using JetNodeConfig = photonjet::JetNodeConfig;
  using JetConfig = photonjet::JetConfig;
  using CentralityConfig = photonjet::CentralityConfig;
  using ProvenanceConfig = photonjet::ProvenanceConfig;

  using SourceRecord = photonjet::SourceRecord;
  using UpstreamRejectedEventRecord = photonjet::UpstreamRejectedEventRecord;
  using EventRecord = photonjet::EventRecord;

  using PhotonTimingRecord = photonjet::PhotonTimingRecord;
  using ShowerShapeRecord = photonjet::ShowerShapeRecord;
  using PhotonCellRecord = photonjet::PhotonCellRecord;
  using IsolationRecord = photonjet::IsolationRecord;
  using IsolationConstituentRecord = photonjet::IsolationConstituentRecord;
  using DominantTruthRecord = photonjet::DominantTruthRecord;
  using PhotonRecord = photonjet::PhotonRecord;
  using NominalPhotonRecord = photonjet::NominalPhotonRecord;
  using NominalReceiptRecord = photonjet::NominalReceiptRecord;

  using JetRecord = photonjet::JetRecord;
  using PhotonJetPairRecord = photonjet::PhotonJetPairRecord;

  using TruthVertexRecord = photonjet::TruthVertexRecord;
  using TruthPhotonRecord = photonjet::TruthPhotonRecord;
  using TruthJetRecord = photonjet::TruthJetRecord;
  using PPG12ClusterRecord = photonjet::PPG12ClusterRecord;
  using EmbeddedDiagnosticRecord = photonjet::EmbeddedDiagnosticRecord;

  using PhotonTruthLinkRecord = photonjet::PhotonTruthLinkRecord;
  using JetTruthLinkRecord = photonjet::JetTruthLinkRecord;

  using WeightComponentRecord = photonjet::WeightComponentRecord;

  using TriggerScalerSnapshotRecord = photonjet::TriggerScalerSnapshotRecord;
  using TriggerRunInfoRecord = photonjet::TriggerRunInfoRecord;

  // ==========================================================================
  // Construction and Fun4All lifecycle
  //
  // PhotonJetTree.cc implements the visible lifecycle:
  //
  //   Init -> InitRun -> process_event -> ResetEvent -> End
  //
  // process_event is the single event-flow spine; the responsibility-specific
  // implementation files supply the individual capture operations it calls.
  // ==========================================================================

  explicit PhotonJetTree(Config config,
                         const std::string& name = "PhotonJetTree");

  ~PhotonJetTree() override;

  PhotonJetTree(const PhotonJetTree&) = delete;
  PhotonJetTree& operator=(const PhotonJetTree&) = delete;

  int Init(PHCompositeNode* topNode) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int ResetEvent(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // ==========================================================================
  // Exposure accounting
  //
  // An event rejected by an upstream reconstruction module never reaches
  // process_event(), but it still belongs to the input exposure. A companion
  // observer registered before reconstruction reports those events through
  // this interface so the original source population remains explicit.
  // ==========================================================================

  /**
   * Persist one event known to have been rejected before this producer ran.
   *
   * This is the only supported path for upstream-rejection accounting; the
   * producer never infers missing encounters from downstream event counts.
   */
  void recordUpstreamRejectedEvent(int run,
                                   std::int64_t physicalEventSequence,
                                   std::int64_t sourceEntry);

  /**
   * Observe one physical input event before reconstruction can reject it.
   *
   * If a new observation arrives while the previous observation is still
   * pending, the previous event never reached process_event() and is recorded
   * as upstream rejected. Implemented in internal/Event.cc.
   */
  void observeUpstreamEvent(std::int64_t physicalEventSequence);

  /**
   * Mark the current product incomplete after an external event-loop failure.
   *
   * End() will then close the file as aborted rather than certifying a partial
   * product as complete.
   */
  void markAborted() { m_aborted = true; }

  /// Read-only source accounting for steering and validation.
  const SourceRecord& sourceRecord() const { return m_source; }

  // Non-owning binding to the registered reconstruction module. Its selected
  // event vertex distinguishes unavailable input from a complete empty photon
  // inventory without a second vertex-selection algorithm in capture.
  void setPhotonBuilder(const photonjet::PhotonClusterBuilder* builder) { m_photonBuilder = builder; }

  // Bind one raw collection to its current-event clustering results.
  // Capture uses exact jet identity and never rematches fractions geometrically.
  using JetFractionReader = std::function<std::array<double, 3>(const Jet*)>;
  void setJetFractionReader(const std::string& node, JetFractionReader reader)
  {
    if (node.empty() || !reader || !m_jetFractionReaders.emplace(node, std::move(reader)).second)
      throw std::invalid_argument("missing or duplicate jet-fraction reader for " + node);
  }

 private:
  // ==========================================================================
  // Production-mode queries
  //
  // Mode is declared configuration. It is never inferred from whether a DST
  // node happens to exist; absence of a required node is a contract failure,
  // not an instruction to switch behavior.
  // ==========================================================================

  bool isPP() const;
  bool isAuAu() const;

  bool isData() const;
  bool isSimulation() const;

  bool isEmbeddedSimulation() const;
  bool isArchivedDoubleInteraction() const;

  bool isPhotonJetSimulation() const;
  bool isInclusiveJetSimulation() const;

  // pp simulation with PPG12 weights keeps the inputs of its PPG12 histogram
  // bookkeeping (PPG12Simulation.cc).
  bool recordsPPG12Simulation() const
  {
    return isSimulation() && isPP() && m_config.weights.ppg12Nominal;
  }

  // ==========================================================================
  // Job and run initialization
  // ==========================================================================

  /**
   * Validate all configuration invariants available before event processing.
   *
   * Contradictory modes, missing provenance, unsupported storage domains, or
   * inconsistent collection definitions fail before a production is allowed
   * to create a seemingly usable event sample.
   */
  void validateConfiguration() const;

  /**
   * Initialize source-level provenance and the compact source identity.
   *
   * Full input hashes remain stored explicitly; the compact identity is the
   * relational key shared by the output tables.
   */
  void initializeSourceRecord();

  /// Open the ROOT product, book its tables, and write file-level provenance.
  void initializeOutput();

  /// Capture state whose lifetime is the current run rather than one event.
  void initializeRun(PHCompositeNode* topNode);

  /// Load configured producer-owned weighting inputs. Implemented in Event.cc.
  void loadWeightInputs();

  /// Load the pinned PPG12 tower mask (pp). Implemented in Objects.cc.
  void loadPPG12TowerMask();

  /**
   * Resolve a final pending upstream observation at end of input.
   *
   * Implemented in internal/Event.cc.
   */
  void finishUpstreamAccounting();

  // ==========================================================================
  // Event context
  //
  // Implemented in internal/Event.cc and called in this order by
  // process_event(). These methods establish the event state to which all
  // reconstructed and truth objects are subsequently attached.
  // ==========================================================================

  void resetEventState();

  void captureEventIdentity(PHCompositeNode* topNode);
  void captureTriggerInformation(PHCompositeNode* topNode);
  void captureVertexInformation(PHCompositeNode* topNode);
  void captureMbdInformation(PHCompositeNode* topNode);
  void captureCentrality(PHCompositeNode* topNode);
  void captureCalorimeterInformation(PHCompositeNode* topNode);
  void captureEventWeights(PHCompositeNode* topNode);

  /// Capture the run-level DATA trigger configuration once per run.
  void captureTriggerRunInformation(PHCompositeNode* topNode);

  /**
   * Record the historical event classification (terminal status, pp
   * simulation sample ownership and the PPG12 audit occurrences; Au+Au
   * embedded stitching, truth-isolation audit and truth-photon misses) after
   * the event context, photons and truth are known. Diagnostic only; no event
   * or object is removed.
   */
  void classifyEvent();
  void classifyAuAuEvent();
  bool embeddedSampleOwnsEvent();
  void recordEmbeddedTruthIsolationAudit();
  void recordEmbeddedTruthMisses();
  bool scaledTriggerFired(const std::string& name) const;

  /*
   * The established producer's PPG12 bookkeeping: the vertex/period audit of
   * a classified event, pp simulation sample ownership and the DATA
   * per-photon yield audit (PPG12.cc); in pp simulation, the histogram fills
   * and occurrences after the weight configuration (PPG12Simulation.cc).
   */
  void recordPPG12VertexAudit(bool preVertexCut);
  bool recordPPG12SampleOwnership(bool beforeEventCuts);
  void recordPPG12PhotonYield();
  void recordPPG12Simulation();

  // ==========================================================================
  // Reconstructed photons
  //
  // Implemented in internal/Objects.cc. This stage stores reconstruction facts,
  // isolation observables and the frozen nominal labels; no photon-ID or
  // isolation working point selects candidates.
  // ==========================================================================

  void capturePhotons(PHCompositeNode* topNode);

  /**
   * Test only the configured photon storage domain.
   *
   * Identification score and isolation decisions deliberately do not enter
   * candidate retention.
   */
  bool photonInCaptureDomain(const PhotonRecord& photon) const;

  void capturePhotonTiming(const RawCluster* cluster,
                           PhotonRecord& photon) const;

  /**
   * Copy builder-produced shower views and retain optional cell witnesses.
   */
  void capturePhotonShowerShapes(const RawCluster* cluster,
                                 TowerInfoContainer* cemcTowers,
                                 PhotonRecord& photon);

  // The isolation cones the builder attached to one candidate.
  std::vector<IsolationRecord> photonIsolation(const RawCluster* cluster,
                                               const PhotonRecord& photon) const;

  void capturePhotonIsolationConstituents(PHCompositeNode* topNode,
                                          const RawCluster* cluster,
                                          const PhotonRecord& photon);

  /**
   * Label every captured candidate with the frozen nominal photon ID and
   * isolation class (implemented in internal/Objects.cc). Labels are stored
   * with the candidate; they never change which candidates are retained.
   */
  void labelNominalPhotons();

  // ==========================================================================
  // Reconstructed jets
  //
  // Implemented in internal/Objects.cc. Collection identity, calibrated/raw jet
  // information, and retained photon-jet relationships are kept explicit.
  // ==========================================================================

  void captureJets(PHCompositeNode* topNode);
  void captureJetCollection(PHCompositeNode* topNode,
                            const JetNodeConfig& collection);

  bool jetInCaptureDomain(const JetRecord& jet) const;

  void buildPhotonJetPairs();

  // ==========================================================================
  // Simulation truth and associations
  //
  // Implemented in internal/Truth.cc. Truth.cc is the single owner of the
  // stored truth population, dominant-primary photon witness, and the final
  // reconstruction-to-truth relations.
  // ==========================================================================

  void captureSimulationTruth(PHCompositeNode* topNode);

  void captureTruthVertices(PHCompositeNode* topNode);
  void captureTruthPhotons(PHCompositeNode* topNode);
  void captureTruthJets(PHCompositeNode* topNode);
  void capturePPG12OwnershipInputs(PHCompositeNode* topNode);

  /**
   * Compute the retained truth-isolation sums and census completeness.
   *
   * A finite cone sum and a complete truth census are independent facts. A bad
   * or unavailable primary may make the census incomplete without replacing a
   * finite measured isolation sum.
   */
  bool calculateTruthPhotonIsolation(PHCompositeNode* topNode,
                                     const PHG4Particle* photon,
                                     double& isolationR03,
                                     double& isolationR04,
                                     bool& complete) const;

  void buildPhotonTruthAssociations(PHCompositeNode* topNode);
  void capturePPG12Clusters(PHCompositeNode* topNode, CaloRawClusterEval& evaluator);
  void buildJetTruthAssociations();

  // ==========================================================================
  // ROOT serialization
  //
  // Implemented in internal/Output.cc. This layer serializes already-defined
  // records and product metadata; it does not introduce new physics quantities.
  // ==========================================================================

  void writeFileMetadata();

  void fillEventOutput();
  void fillRejectedOutput(const UpstreamRejectedEventRecord& record);

  void writeCompletionMetadata();

  /// Write the product and release all ROOT output handles.
  void closeOutput();

  // ==========================================================================
  // Immutable job configuration
  // ==========================================================================

  Config m_config;
  const photonjet::PhotonClusterBuilder* m_photonBuilder = nullptr;  // owned by Fun4AllServer

  // The frozen nominal photon-ID model, loaded once per job.
  std::shared_ptr<const photonjet::production::PinnedModel> m_nominalModel;
  std::map<std::string, JetFractionReader> m_jetFractionReaders;

  // ==========================================================================
  // ROOT product handles
  //
  // One bounded direct writer. No intermediate object TTrees are created.
  // ==========================================================================

  std::shared_ptr<photonjet::DirectEventWriter> m_directOutput;

  // ==========================================================================
  // Source- and run-lifetime state
  //
  // These values persist across ResetEvent(). They describe the source as a
  // whole, upstream exposure accounting, and run-level trigger information.
  // ==========================================================================

  SourceRecord m_source;

  std::vector<UpstreamRejectedEventRecord> m_upstreamRejectedEvents;

  // Most recent input observation not yet consumed by process_event().
  bool m_upstreamPending = false;
  int m_upstreamPendingRun = 0;
  std::int64_t m_upstreamPendingSequence = -1;
  std::uint64_t m_upstreamObservedEvents = 0;

  std::vector<TriggerScalerSnapshotRecord> m_triggerScalerSnapshots;
  std::vector<TriggerRunInfoRecord> m_triggerRunInfo;

  // Source totals of the nominal labels, written with the file's receipts.
  NominalReceiptRecord m_nominalReceipt;

  // Producer-owned simulated-vertex reweighting histogram.
  TH1* m_vertexReweight = nullptr;

  // PPG12 tower mask for the photon-yield source arbitration (pp only).
  TH2* m_ppg12TowerMask = nullptr;

  // ==========================================================================
  // Current-event state
  //
  // resetEventState() clears these containers for every producer encounter.
  // Source- and run-lifetime state above is unaffected.
  // ==========================================================================

  EventRecord m_event;

  std::vector<PhotonRecord> m_photons;
  std::vector<PhotonCellRecord> m_photonCells;
  std::vector<IsolationRecord> m_isolationRecords;
  std::vector<IsolationConstituentRecord> m_isolationConstituents;

  std::vector<JetRecord> m_jets;
  std::vector<PhotonJetPairRecord> m_photonJetPairs;

  std::vector<TruthVertexRecord> m_truthVertices;
  std::vector<TruthPhotonRecord> m_truthPhotons;
  std::vector<TruthJetRecord> m_truthJets;

  // The identity association, and the established producer's distance links.
  std::vector<PhotonTruthLinkRecord> m_photonTruthLinks;
  std::vector<PhotonTruthLinkRecord> m_photonDistanceLinks;
  std::vector<JetTruthLinkRecord> m_jetTruthLinks;

  std::vector<WeightComponentRecord> m_weightComponents;

  // Inputs of the PPG12 simulation bookkeeping, not stored: every photon
  // candidate in container order, the builder's input clusters, the unsplit
  // response clusters when that node exists, the radius codes of the truth
  // jet containers read, and whether the cluster truth evaluator had its nodes.
  std::vector<PhotonRecord> m_ppg12Photons;
  std::vector<PPG12ClusterRecord> m_ppg12Clusters;
  std::vector<PPG12ClusterRecord> m_ppg12ResponseClusters;
  bool m_ppg12ResponseClustersPresent = false;
  std::vector<int> m_ppg12TruthJetRadiiFound;
  bool m_ppg12EvaluatorReady = false;

  // ==========================================================================
  // Lifecycle guards
  //
  // These flags encode producer progress independently of ROOT pointer state.
  // In particular, an opened file is not equivalent to a completed product.
  // ==========================================================================

  bool m_outputReady = false;
  bool m_runInitialized = false;
  bool m_finalized = false;
  bool m_aborted = false;
};

#endif  // PHOTONJETTREE_H
