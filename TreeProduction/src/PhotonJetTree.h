#ifndef PHOTONJETTREE_H
#define PHOTONJETTREE_H

/**
 * @file PhotonJetTree.h
 * @brief Public interface and state of the isolated-photon + jet tree producer.
 *
 * PhotonJetTree is the single Fun4All producer used for pp and AuAu, DATA and simulation. Collision-system and sample differences are expressed through
 * typed configuration rather than through separate producer classes.
 *
 * The producer owns the conversion of reconstructed event state into canonical
 * analysis records:
 *
 *   source / event identity and exposure accounting
 *     -> trigger and scaler state
 *     -> vertex, MBD, centrality, and calorimeter state
 *     -> reconstructed photons, shower information, and isolation
 *     -> reconstructed jets and photon-jet relationships
 *     -> simulation truth and reco-truth associations
 *     -> producer-owned event weights
 *     -> ROOT serialization and provenance
 *
 * TreeProduction deliberately stops at reconstructed and truth-level facts.
 * Detector reconstruction is registered by the production layer; photon-ID models, working points, analysis selections, histogramming, corrections, unfolding, normalization, and plotting are downstream responsibilities.
 *
 * Candidates are therefore retained independently of photon-ID and isolation working points. The trees preserve the reconstruction primitives needed to train or evaluate later classifiers without rereading the DST.
 *
 * The implementation is split by responsibility while remaining one class:
 *
 *   PhotonJetTree.cc          lifecycle and explicit event-flow spine
 *   internal/Event.cc        event context, exposure, triggers, and weights
 *   internal/Photons.cc      photons, shower information, and isolation
 *   internal/Jets.cc         reconstructed jets and photon-jet relations
 *   internal/Truth.cc        simulation truth and reco-truth associations
 *   internal/Output.cc       ROOT serialization, provenance, and completion
 *
 * All persisted configuration and record types are defined once in
 * internal/Types.h.
 */

#include "internal/Types.h"

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>
#include <vector>

class PHCompositeNode;

class TFile;
class TH1;
class TTree;

class RawCluster;
class TowerInfoContainer;
class PhotonClusterBuilder;

class PHG4Particle;

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

  using JetRecord = photonjet::JetRecord;
  using PhotonJetPairRecord = photonjet::PhotonJetPairRecord;

  using TruthVertexRecord = photonjet::TruthVertexRecord;
  using TruthPhotonRecord = photonjet::TruthPhotonRecord;
  using TruthJetRecord = photonjet::TruthJetRecord;

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
  void observeUpstreamEvent(int run, std::int64_t physicalEventSequence);

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
  void setPhotonBuilder(const PhotonClusterBuilder* builder) { m_photonBuilder = builder; }

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

  // ==========================================================================
  // Reconstructed photons
  //
  // Implemented in internal/Photons.cc. This stage stores reconstruction facts
  // and isolation observables; it does not apply photon-ID or isolation working
  // points.
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

  void capturePhotonIsolation(const RawCluster* cluster,
                              PhotonRecord& photon);

  void capturePhotonIsolationConstituents(PHCompositeNode* topNode,
                                          const RawCluster* cluster,
                                          const PhotonRecord& photon);

  // ==========================================================================
  // Reconstructed jets
  //
  // Implemented in internal/Jets.cc. Collection identity, calibrated/raw jet
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
  void buildJetTruthAssociations();

  // ==========================================================================
  // ROOT serialization
  //
  // Implemented in internal/Output.cc. This layer serializes already-defined
  // records and product metadata; it does not introduce new physics quantities.
  // ==========================================================================

  void bookTrees();
  void writeFileMetadata();

  void fillEventOutput();
  void fillObjectOutput();

  void writeCompletionMetadata();

  /// Write the product and release all ROOT output handles.
  void closeOutput();

  // ==========================================================================
  // Immutable job configuration
  // ==========================================================================

  Config m_config;
  const PhotonClusterBuilder* m_photonBuilder = nullptr;  // owned by Fun4AllServer

  // ==========================================================================
  // ROOT product handles
  //
  // One handle per persisted table. The corresponding row buffers are declared
  // separately below.
  // ==========================================================================

  TFile* m_outputFile = nullptr;

  TTree* m_sourceTree = nullptr;
  TTree* m_upstreamRejectedTree = nullptr;
  TTree* m_eventTree = nullptr;

  TTree* m_photonTree = nullptr;
  TTree* m_showerShapeTree = nullptr;
  TTree* m_photonCellTree = nullptr;
  TTree* m_isolationTree = nullptr;
  TTree* m_isolationConstituentTree = nullptr;

  TTree* m_jetTree = nullptr;
  TTree* m_photonJetPairTree = nullptr;

  TTree* m_truthVertexTree = nullptr;
  TTree* m_truthPhotonTree = nullptr;
  TTree* m_truthJetTree = nullptr;

  TTree* m_photonTruthLinkTree = nullptr;
  TTree* m_jetTruthLinkTree = nullptr;

  TTree* m_weightTree = nullptr;

  TTree* m_triggerScalerTree = nullptr;
  TTree* m_triggerRunInfoTree = nullptr;

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

  // Producer-owned simulated-vertex reweighting histogram.
  TH1* m_vertexReweight = nullptr;

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

  std::vector<PhotonTruthLinkRecord> m_photonTruthLinks;
  std::vector<JetTruthLinkRecord> m_jetTruthLinks;

  std::vector<WeightComponentRecord> m_weightComponents;

  // ==========================================================================
  // ROOT row buffers
  //
  // ROOT branches hold addresses into these records. Output.cc copies the row
  // being serialized into the appropriate buffer before each TTree::Fill().
  // ==========================================================================

  SourceRecord m_sourceRow;
  UpstreamRejectedEventRecord m_upstreamRejectedRow;
  EventRecord m_eventRow;

  PhotonRecord m_photonRow;
  ShowerShapeRecord m_showerShapeRow;
  PhotonCellRecord m_photonCellRow;
  IsolationRecord m_isolationRow;
  IsolationConstituentRecord m_isolationConstituentRow;

  JetRecord m_jetRow;
  PhotonJetPairRecord m_photonJetPairRow;

  TruthVertexRecord m_truthVertexRow;
  TruthPhotonRecord m_truthPhotonRow;
  TruthJetRecord m_truthJetRow;

  PhotonTruthLinkRecord m_photonTruthLinkRow;
  JetTruthLinkRecord m_jetTruthLinkRow;

  WeightComponentRecord m_weightRow;

  TriggerScalerSnapshotRecord m_triggerScalerRow;
  TriggerRunInfoRecord m_triggerRunInfoRow;

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
