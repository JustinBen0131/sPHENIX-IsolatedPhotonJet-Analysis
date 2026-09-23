#ifndef PHOTONJETTREE_H
#define PHOTONJETTREE_H

//
// PhotonJetTree
//
// One Fun4All module that reads reconstructed sPHENIX events and writes the
// isolated-photon + jet analysis trees. One producer serves p+p and Au+Au,
// data and simulation; the differences are explicit typed configuration, not
// separate classes.
//
//
// What this module owns
//
//
//   source and event identity, and the exposure accounting that goes with it
//   trigger decisions and scaler bookkeeping
//   reconstructed vertex, minimum-bias detector state, centrality
//   calorimeter event state
//   photon candidates, their timing, shower shapes, cells and isolation
//   reconstructed jets and their photon-jet relationships
//   simulation truth, and the association of truth to reconstruction
//   producer event weights
//   the ROOT output and its provenance
//
//
// What this module does not own
//
//
//   detector reconstruction and calibration      the production layer
//                                                registers those modules
//   photon identification models                 trained and applied
//                                                downstream, on these trees
//   working points, isolation and identification selections
//   background subtraction, purity, response matrices, unfolding
//   histograms, normalisation and plotting
//
// The module is model-independent by construction. Candidates are admitted on
// kinematics alone, and the trees carry the complete reconstruction primitives
// a classifier needs, so a new model can be trained and applied without
// reading a DST again.
//
//
// Event flow
//
//
//   identity -> trigger -> vertex -> minimum bias -> centrality -> calorimeter
//            -> photons -> jets -> truth -> associations -> weights -> output
//
// process_event follows exactly that order and nothing hides behind a callback
// table. The typed records it fills are defined once in internal/Types.h.
//
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

class PHG4Particle;

class PhotonJetTree : public SubsysReco
{
 public:
  //
  // The canonical vocabulary, defined in internal/Types.h.
  //
  // These are names for existing types, not second definitions. They let the
  // implementation files spell a record the way a reader expects to see it on
  // this class.
  //
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

  //
  // Construction and Fun4All lifecycle
  //
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

  //
  // Exposure accounting for events this module never saw
  //
  // When a reconstruction module aborts an event, Fun4All skips every module
  // after it, so this producer is never called. Those events are still
  // exposure. A companion observer registered ahead of reconstruction reports
  // them here, and they are written to their own table and counted in the
  // source record.
  //
  // Calling this is the only supported way to account for them. The producer
  // will not invent the count.
  //
  void recordUpstreamRejectedEvent(int run,
                                   std::int64_t physicalEventSequence,
                                   std::int64_t sourceEntry);

  // Called by the exposure observer the production layer registers ahead of
  // reconstruction, once per input event, before any module can abort it. If
  // a second observation arrives while the previous one has not been consumed
  // by process_event, the previous event never reached this producer and is
  // recorded as upstream rejected. Implemented in internal/Event.cc.
  void observeUpstreamEvent(int run, std::int64_t physicalEventSequence);

  // Called by the steering macro when the Fun4All event loop stopped early.
  // End() then closes the file with completion_status "aborted" instead of
  // "complete", so a partial file can never pass for a finished one.
  void markAborted() { m_aborted = true; }

  // Read-only view of the running source accounting, for validation code.
  const SourceRecord& sourceRecord() const { return m_source; }

 private:
  //
  // Production mode
  //
  // Mode is configuration, never inferred from whether a node happens to be
  // present. A missing node in a mode that requires it is an error, not a
  // signal to switch behaviour.
  //
  bool isPP() const;
  bool isAuAu() const;

  bool isData() const;
  bool isSimulation() const;

  bool isEmbeddedSimulation() const;
  bool isArchivedDoubleInteraction() const;

  bool isPhotonJetSimulation() const;
  bool isInclusiveJetSimulation() const;

  //
  // Setup
  //
  // Reject an internally inconsistent job before any file is opened. Throws
  // std::runtime_error describing the first inconsistency found.
  void validateConfiguration() const;

  // Fill the source record and its identity from the configuration. The
  // identity needs a pinned input digest and does not depend on any event.
  void initializeSourceRecord();

  // Open the output file, book every tree, and write the file-level
  // provenance. Implemented in internal/Output.cc.
  void initializeOutput();

  // Capture run-dependent state, such as the trigger configuration table.
  void initializeRun(PHCompositeNode* topNode);

  // Open the vertex-reweighting histogram when the configuration asks for
  // it. Implemented in internal/Event.cc.
  void loadWeightInputs();

  // Resolve a pending upstream observation that process_event never
  // consumed, at the end of the event stream. Implemented in internal/Event.cc.
  void finishUpstreamAccounting();

  //
  // Per-event capture, in event order
  //
  // Implemented in internal/Event.cc.
  //
  void resetEventState();

  void captureEventIdentity(PHCompositeNode* topNode);
  void captureTriggerInformation(PHCompositeNode* topNode);
  void captureVertexInformation(PHCompositeNode* topNode);
  void captureMbdInformation(PHCompositeNode* topNode);
  void captureCentrality(PHCompositeNode* topNode);
  void captureCalorimeterInformation(PHCompositeNode* topNode);
  void captureEventWeights(PHCompositeNode* topNode);

  // Run-level trigger configuration table, read once per run in data.
  void captureTriggerRunInformation(PHCompositeNode* topNode);

  //
  // Reconstructed photons
  //
  // Implemented in internal/Photons.cc.
  //
  void capturePhotons(PHCompositeNode* topNode);

  // Storage domain only: finite kinematics inside the configured transverse
  // energy and pseudorapidity window. No identification and no isolation.
  bool photonInCaptureDomain(const PhotonRecord& photon) const;

  void capturePhotonTiming(const RawCluster* cluster,
                           TowerInfoContainer* cemcTowers,
                           PhotonRecord& photon) const;

  // Evaluate every configured shower-shape definition from the calibrated
  // tower grid around the cluster, and record the cells that entered it.
  void capturePhotonShowerShapes(const RawCluster* cluster,
                                 TowerInfoContainer* cemcTowers,
                                 PhotonRecord& photon);

  void capturePhotonIsolation(const RawCluster* cluster,
                              PhotonRecord& photon);

  void capturePhotonIsolationConstituents(PHCompositeNode* topNode,
                                          const RawCluster* cluster,
                                          const PhotonRecord& photon);

  //
  // Reconstructed jets
  //
  // Implemented in internal/Jets.cc.
  //
  void captureJets(PHCompositeNode* topNode);
  void captureJetCollection(PHCompositeNode* topNode,
                            const JetNodeConfig& collection);

  bool jetInCaptureDomain(const JetRecord& jet) const;

  void buildPhotonJetPairs();

  //
  // Simulation truth and associations
  //
  // Implemented in internal/Truth.cc. Truth.cc is the single owner of the
  // dominant-primary witness as well as of the final relations.
  //
  void captureSimulationTruth(PHCompositeNode* topNode);

  void captureTruthVertices(PHCompositeNode* topNode);
  void captureTruthPhotons(PHCompositeNode* topNode);
  void captureTruthJets(PHCompositeNode* topNode);

  // Finite isolation and complete census are separate historical facts. A
  // skipped bad primary marks completeness false without changing a finite
  // stored sum or the native signal predicate.
  bool calculateTruthPhotonIsolation(PHCompositeNode* topNode,
                                     const PHG4Particle* photon,
                                     double& isolationR03,
                                     double& isolationR04,
                                     bool& complete) const;

  void buildPhotonTruthAssociations(PHCompositeNode* topNode);
  void buildJetTruthAssociations();

  //
  // ROOT output
  //
  // Implemented in internal/Output.cc. No scientific quantity is computed for
  // the first time there.
  //
  void bookTrees();
  void writeFileMetadata();

  void fillEventOutput();
  void fillObjectOutput();

  void writeCompletionMetadata();

  // Write and close the file, releasing the tree handles.
  void closeOutput();

  //
  // Configuration
  //
  Config m_config;

  //
  // Output file and tables
  //
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

  //
  // State that lives as long as the source
  //
  SourceRecord m_source;

  std::vector<UpstreamRejectedEventRecord> m_upstreamRejectedEvents;

  // The most recent upstream observation not yet consumed by process_event.
  bool m_upstreamPending = false;
  int m_upstreamPendingRun = 0;
  std::int64_t m_upstreamPendingSequence = -1;
  std::uint64_t m_upstreamObservedEvents = 0;

  std::vector<TriggerScalerSnapshotRecord> m_triggerScalerSnapshots;
  std::vector<TriggerRunInfoRecord> m_triggerRunInfo;

  // Simulated-vertex reweighting histogram, owned by this module.
  TH1* m_vertexReweight = nullptr;

  //
  // State of the event being processed
  //
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

  //
  // Row buffers
  //
  // ROOT branches address these; a fill copies the record being written into
  // the matching buffer first.
  //
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

  //
  // Lifecycle guards
  //
  bool m_outputReady = false;
  bool m_runInitialized = false;
  bool m_finalized = false;
  bool m_aborted = false;
};

#endif  // PHOTONJETTREE_H
