#ifndef PHOTONJETTREE_INTERNAL_PRODUCTION_H
#define PHOTONJETTREE_INTERNAL_PRODUCTION_H

/**
 * @file Production.h
 * Shared Job, Plan, conditions and setup interface.
 * The frontend supplies source identity; Production.cc resolves the preset.
 * Setup.cc binds conditions, reconstruction and physical input streams.
 * Record definitions are owned by Types.h.
 */

#include "Types.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

class Fun4AllServer;
class PhotonJetTree;

namespace photonjet
{
namespace production
{

// Read native rule coverage before constructing an output or running events.
// The resolved schema is authoritative; only shared nominal augmentation may
// remain missing when ready is true. This is implementation readiness only.
std::string OutputMappingStatusJson(const Config& config);
// Atomic progress from the original controlling process, in outputFile +
// ".progress.json". Negative processed arguments reuse its last observation.
// A finite source budget is mandatory; invalid/regressing counts throw.
void WriteProgress(const Config& config,const std::string& phase,
    std::int64_t rowProcessed,std::int64_t skipProcessed,
    const std::string& error="",bool force=false);

// ============================================================================
// Physical input layouts
//
// InputLayout describes how the physical files belonging to one source are
// presented to Fun4All. It is an operational property of the source, not a
// physics label written by PhotonJetTree.
// ============================================================================

enum class InputLayout
{
  PairedData,                 // calorimeter file plus an auxiliary file, synchronised
  ReconstructedSimulation,    // calorimeter, Geant4, truth jets, global, MBD lanes
  ArchivedDoubleInteraction   // Geant4 hits plus prebuilt truth jets only
};


// ============================================================================
// Truth-jet source
//
// Simulation may consume truth jets shipped with the source or construct them
// during this production pass. Auto is resolved only after the physical input
// inventory is known.
// ============================================================================

enum class TruthJetMode
{
  None,   // data
  Auto,   // DST when the input carries them, otherwise Build
  DST,    // read the truth-jet collections shipped with the input
  Build   // cluster them from generator particles during this pass
};


// ============================================================================
// Per-source production job
//
// Job contains values that vary for every invocation and therefore do not
// belong in the shared tree_production.yaml configuration.
//
// The production frontend constructs this object after resolving one exact
// source from the reviewed population manifest. LoadPlan combines it with one
// named YAML profile to produce the complete Plan.
// ============================================================================

struct Job
{
  std::string profile;

  std::string inputList;
  std::string outputFile;

  int run = 0;
  int segment = 0;

  std::int64_t firstEntry = 0;
  std::int64_t nEvents = -1;

  std::int64_t sourceFileOrdinal = -1;

  std::string dataset;
  std::string sample;
  std::string period;
  std::string siDiRole;

  std::string inputUriSha256;
  std::string inputFileSha256;
  std::string sourceManifestSha256;

  // Hash of the actual steering macro, supplied by that macro itself.
  std::string macroSha256;

  int verbosity = 0;
};


// ============================================================================
// One Fun4All input stream
//
// A physical source can require several logical streams that advance together.
// The stream name identifies its reconstruction role; files contains the
// already-resolved physical input.
//
// synchronized=true uses ordinary Fun4All synchronization. NoSync is an
// explicit production decision for modes whose inputs must advance
// independently, notably embedded and archived sources.
// ============================================================================

struct InputStream
{
  std::string name;
  std::vector<std::string> files;

  // Ordinary Fun4All synchronisation is the default. NoSync is an explicit
  // production decision for embedded and archived inputs.
  bool synchronized = true;

  bool required = true;
};


// ============================================================================
// Run-dependent conditions
//
// Conditions contains the external calibration state required to instantiate
// one Plan. Production.cc resolves the intended identities; ConfigureConditions
// binds CDB/local payload locations and verifies pinned digests before the event
// loop begins.
//
// Paths in this structure refer to the exact resolved artifacts used by the
// production. Missing required conditions are never silently replaced.
// ============================================================================

struct Conditions
{
  std::string cdbGlobalTag;
  std::uint64_t cdbTimestamp = 0;


  // --------------------------------------------------------------------------
  // Au+Au DATA centrality
  //
  // Local mode binds a coherent divisions/run-scale/vertex-scale triplet to a
  // specific run. Own-run payloads are the default. The only permitted
  // cross-run fallback is separately qualified and requires pinned evidence.
  // --------------------------------------------------------------------------

  bool centralityLocalPayloads = false;

  std::string centralityCalibrationSet;
  std::string centralityManifestSha256;

  int centralityPayloadRun = 0;

  std::string centralityDivisions;
  std::string centralityDivisionsSha256;

  std::string centralityRunScale;
  std::string centralityRunScaleSha256;

  std::string centralityVertexScale;
  std::string centralityVertexScaleSha256;

  std::string centralityFallbackEvidenceSha256;


  // --------------------------------------------------------------------------
  // CEMC tower status, Au+Au DATA
  //
  // The paired calorimeter input carries a historical tower-status state.
  // The guard verifies that state against the bad-tower map used by the
  // original reconstruction.
  //
  // If that historical map was not persisted, recovery is allowed only through
  // an explicitly authorized, run-bound, digest-pinned payload.
  // --------------------------------------------------------------------------

  bool cemcStatusGuard = false;

  std::string cemcOriginalCalorimeterInput;

  bool cemcStatusRecoveryAuthorized = false;

  std::string cemcStatusRecoveryPayload;
  std::string cemcStatusRecoveryPayloadSha256;


  // --------------------------------------------------------------------------
  // CEMC zero-suppression cross calibration
  //
  // DATA only. The resolved payload is pinned and validated before the CEMC
  // calibration module executes.
  // --------------------------------------------------------------------------

  bool cemcZeroSuppressionCrossCalibration = false;

  std::string cemcZeroSuppressionPayload;


  // --------------------------------------------------------------------------
  // Jet energy scale
  //
  // JES is applied upstream exactly once. PhotonJetTree records the resulting
  // corrected and raw jet quantities together with this payload identity.
  // --------------------------------------------------------------------------

  bool applyJetEnergyScale = true;

  std::string jetEnergyScaleCdbKey = "JES_Calib_Default";
  std::string jetEnergyScalePayload;
  std::string jetEnergyScalePayloadSha256;

  bool jetEnergyScaleUsesEmFraction = false;


  // --------------------------------------------------------------------------
  // Simulated-vertex reweighting
  //
  // This is producer-owned because it is intrinsic to the simulated event
  // record. Cross-section, stitching, centrality, SI/DI, and luminosity weights
  // remain downstream.
  // --------------------------------------------------------------------------

  bool applyVertexReweight = false;

  std::string vertexReweightFile;
  std::string vertexReweightFileSha256;
  std::string vertexReweightHistogram;
};


// ============================================================================
// Reconstruction steering
//
// Reconstruction is the resolved upstream reconstruction graph for one profile:
// which stages run and which parameters they receive.
//
// These fields configure Fun4All registration. They are not event observables
// and are not a substitute for the measurements PhotonJetTree records.
// ============================================================================

struct Reconstruction
{
  InputLayout inputLayout = InputLayout::PairedData;
  TruthJetMode truthJetMode = TruthJetMode::None;


  // --------------------------------------------------------------------------
  // Production mode
  //
  // These flags are resolved once from the selected profile and then drive all
  // later reconstruction decisions.
  // --------------------------------------------------------------------------

  bool isPP = false;
  bool isAuAu = false;

  bool isData = false;
  bool isSimulation = false;

  bool isEmbedded = false;
  bool isArchivedDoubleInteraction = false;


  // --------------------------------------------------------------------------
  // Reconstruction stages
  // --------------------------------------------------------------------------

  // Minimum-bias detector and global vertex reconstruction. Off for embedded
  // simulation, which consumes the products of its minimum-bias event.
  bool reconstructEventGeometry = false;

  // Tower calibration, status and clustering from raw or simulated towers.
  // Off for reconstructed p+p simulation, whose input is already calibrated.
  bool reconstructCalorimeter = false;

  bool runCalorimeterStatusSkimmer = false;  // data: drop incomplete events
  bool setTowerStatus = false;               // legacy status setters
  bool reconstructZdc = false;               // data: needed by minimum bias

  bool reconstructCentrality = false;
  bool runMinimumBiasClassifier = false;

  bool reconstructAuAuBackground = false;

  bool reconstructJets = true;
  bool reconstructPhotons = true;


  // --------------------------------------------------------------------------
  // Calorimeter reconstruction
  // --------------------------------------------------------------------------

  double clusterThresholdGeV = 0.070;
  std::string clusterProfile;

  // Simulation applies a hot-tower map from the conditions database.
  std::string simulationHotTowerMapKey = "CEMC_hotTowers_status";


  // --------------------------------------------------------------------------
  // Photon builder
  //
  // These are reconstruction-admission settings. They are deliberately aligned
  // with the broad PhotonJetTree storage domain so the builder does not remove
  // candidates that downstream analysis may later need.
  // --------------------------------------------------------------------------

  double photonMinEtGeV = 5.0;

  bool photonUseVertexCut = false;
  double photonVertexCutCm = 60.0;

  double showerShapeTowerMinEnergyGeV = 0.070;

  bool useTopoclusterIsolation = false;
  std::string topoclusterNode = "TOPOCLUSTER_ALLCALO";

  bool useSubtractedIsolation = false;


  // --------------------------------------------------------------------------
  // Au+Au underlying-event subtraction
  // --------------------------------------------------------------------------

  double retowerFractionCut = 0.5;

  // 0 no flow, 1 psi2 from the calorimeter, 2 from HIJING, 3 from the sEPD.
  int flowModulation = 0;

  int firstPassSeedJetD = 3;
  double secondPassSeedPtGeV = 7.0;


  // --------------------------------------------------------------------------
  // Archived p+p double-interaction reconstruction
  //
  // The archived Geant4-only donor production consumed a fixed random-seed
  // sequence in order. Exact reproduction therefore requires FIFO replay of
  // that sequence and forbids an independent RecoConsts seed from silently
  // perturbing the historical random-number stream.
  // --------------------------------------------------------------------------

  std::vector<unsigned int> historicalSeeds;

  int expectedPedestalSequence = -1;

  bool forbidRecoConstsRandomSeed = true;
};


// ============================================================================
// Fully resolved production plan
//
// Plan is the complete contract handed to the steering macro after LoadPlan and
// Validate. Once constructed, the macro does not need to understand YAML,
// profile inheritance, node names, input-list structure, or payload policy.
//
// The exact configuration bytes and digest are carried alongside the typed
// state so they can be embedded in output provenance.
// ============================================================================

struct Plan
{
  Job job;

  photonjet::Config tree;

  std::vector<InputStream> inputs;

  Conditions conditions;
  Reconstruction reconstruction;

  // The exact configuration bytes used to build this plan.
  std::string configurationText;
  std::string configurationSha256;
};


// ============================================================================
// Shared implementation details
//
// Only helpers required across the Production translation units are declared
// here. File-local implementation details should remain private to their
// owning .cc file so this header stays a compact steering contract.
// ============================================================================

namespace detail
{

// Historical conditions run used by the archived double-interaction source.
constexpr int kArchivedDoubleInteractionRun = 28;

// Run numbers above this threshold are treated as real DATA-run coordinates by
// the current production contract.
constexpr int kDataRunThreshold = 1000;

// Number of CEMC channels expected by status/calibration validation.
constexpr std::size_t kCemcChannelCount = 24576;


/**
 * Throw a consistently prefixed production error.
 */
[[noreturn]] void fail(const std::string& message);


/**
 * Remove leading and trailing ASCII whitespace.
 */
std::string trim(const std::string& value);


/**
 * Convert a configuration token to lowercase for controlled parsing.
 */
std::string lower(std::string value);


/**
 * Return true only for a syntactically valid lowercase/hex SHA256 digest.
 */
bool validSha256(const std::string& digest);


/**
 * Compute the SHA256 digest of one local file.
 */
std::string fileSha256(const std::string& path);


/**
 * Interpret the physical input list for one resolved Job and reconstruction
 * mode.
 *
 * Implemented by Setup.cc.
 */
std::vector<InputStream> loadInputs(const Job& job,
                                    const Reconstruction& reconstruction);


/**
 * Test whether a named resolved input stream exists in the Plan.
 */
bool hasInputStream(const Plan& plan,
                    const std::string& name);


/**
 * Convert an upstream status/return-code predicate into a production failure.
 */
void requireStatus(bool condition,
                   const std::string& reason);

}  // namespace detail


// ============================================================================
// Procedural production interface
//
// Fun4All_PhotonJetTree.C uses only these functions. Their order mirrors the
// lifetime of one production job:
//
//   LoadPlan
//       resolve Job + YAML + physical inputs into typed steering state
//
//   ConfigureConditions
//       bind and verify run-dependent external payloads
//
//   Print
//       expose the major resolved choices before the event loop
//
//   RegisterReconstruction
//       register observer, reconstruction graph, and PhotonJetTree
//
//   RegisterInputs
//       attach the physical Fun4All input streams
//
// Validate is public because LoadPlan uses it and because resolved plans can be
// checked independently before execution.
// ============================================================================

/**
 * Parse the shared configuration, resolve one named profile and the physical
 * input list, construct the typed PhotonJetTree/reconstruction/conditions
 * state, and validate the resulting Plan.
 */
Plan LoadPlan(const std::string& yamlFile,
              const Job& job);


/**
 * Reject an internally inconsistent or scientifically unbound Plan before
 * production begins.
 */
void Validate(const Plan& plan);


/**
 * Bind the run context, resolve CDB/local conditions payloads, verify every
 * required digest, and write the resolved payload paths back into the Plan.
 */
void ConfigureConditions(Plan& plan);


/**
 * Register the upstream exposure observer, the reconstruction modules required
 * by the resolved profile, and finally PhotonJetTree.
 */
void RegisterReconstruction(Fun4AllServer* server,
                            const Plan& plan,
                            PhotonJetTree* producer);


/**
 * Register the resolved physical input streams and, when required by
 * reconstruction in this pass, the geometry run node.
 */
void RegisterInputs(Fun4AllServer* server,
                    const Plan& plan);


/**
 * Print a concise human-readable summary of the resolved production job.
 */
void Print(const Plan& plan);

}  // namespace production
}  // namespace photonjet

#endif  // PHOTONJETTREE_INTERNAL_PRODUCTION_H
