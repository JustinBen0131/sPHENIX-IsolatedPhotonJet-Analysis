#ifndef PHOTONJETTREE_INTERNAL_PRODUCTION_H
#define PHOTONJETTREE_INTERNAL_PRODUCTION_H

//
// Production.h
//
// How the reconstructed event that PhotonJetTree reads is created.
//
// The producer captures what it finds on the node tree. Everything that puts
// it there, from calorimeter calibration to jet energy scale, is registered by
// this layer from one resolved Plan. The macro sees only the six functions at
// the bottom of this file and never touches YAML, node names or module
// parameters itself.
//
// Nothing here computes an event observable.
//
#include "Types.h"

#include <cstdint>
#include <string>
#include <vector>

class Fun4AllServer;
class PhotonJetTree;

namespace photonjet
{
namespace production
{

// How the physical input list is laid out and synchronised.
enum class InputLayout
{
  PairedData,                 // calorimeter file plus an auxiliary file, synchronised
  ReconstructedSimulation,    // calorimeter, Geant4, truth jets, global, MBD lanes
  ArchivedDoubleInteraction   // Geant4 hits plus prebuilt truth jets only
};

// Where truth jets come from in simulation.
enum class TruthJetMode
{
  None,   // data
  Auto,   // DST when the input carries them, otherwise Build
  DST,    // read the truth-jet collections shipped with the input
  Build   // cluster them from generator particles during this pass
};

// One production invocation. These values change from job to job and are
// never written into the shared configuration file.
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

  int verbosity = 0;
};

// One Fun4All input stream.
struct InputStream
{
  std::string name;
  std::vector<std::string> files;

  // Ordinary Fun4All synchronisation is the default. NoSync is an explicit
  // production decision for embedded and archived inputs.
  bool synchronized = true;

  bool required = true;
};

// Run-dependent conditions and the payloads they resolve to. Every path here
// is verified against its pinned digest before the event loop starts.
struct Conditions
{
  std::string cdbGlobalTag;
  std::uint64_t cdbTimestamp = 0;

  // ---- Au+Au data centrality -------------------------------------------
  //
  // Local mode binds the run's own frozen divisions, run scale and vertex
  // scale. A missing own-run payload never falls back to a neighbouring run;
  // the one qualified fallback run needs pinned evidence.

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

  // ---- CEMC tower status, Au+Au data -----------------------------------
  //
  // Paired Au+Au data carries its tower status through the calibration. The
  // guard validates that status against the bad-tower map the original file
  // was reconstructed with. When the original file has no saved map, an
  // explicitly authorised, run-bound recovery payload restores the hot flags.

  bool cemcStatusGuard = false;

  std::string cemcOriginalCalorimeterInput;

  bool cemcStatusRecoveryAuthorized = false;
  std::string cemcStatusRecoveryPayload;
  std::string cemcStatusRecoveryPayloadSha256;

  // CEMC zero-suppression cross calibration, data only. The resolved payload
  // is pinned and its ratios validated before calibration runs.
  bool cemcZeroSuppressionCrossCalibration = false;
  std::string cemcZeroSuppressionPayload;

  // ---- Jet energy scale ------------------------------------------------

  bool applyJetEnergyScale = true;

  std::string jetEnergyScaleCdbKey = "JES_Calib_Default";
  std::string jetEnergyScalePayload;
  std::string jetEnergyScalePayloadSha256;

  bool jetEnergyScaleUsesEmFraction = false;

  // ---- Simulated-vertex reweighting ------------------------------------

  bool applyVertexReweight = false;

  std::string vertexReweightFile;
  std::string vertexReweightFileSha256;
  std::string vertexReweightHistogram;
};

// The reconstruction graph for one profile: which modules run and with what
// parameters. Steering state only.
struct Reconstruction
{
  InputLayout inputLayout = InputLayout::PairedData;
  TruthJetMode truthJetMode = TruthJetMode::None;

  // ---- Mode ------------------------------------------------------------

  bool isPP = false;
  bool isAuAu = false;

  bool isData = false;
  bool isSimulation = false;

  bool isEmbedded = false;
  bool isArchivedDoubleInteraction = false;

  // ---- Stages ----------------------------------------------------------

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

  // ---- Calorimeter -----------------------------------------------------

  double clusterThresholdGeV = 0.070;
  std::string clusterProfile;

  // Simulation applies a hot-tower map from the conditions database.
  std::string simulationHotTowerMapKey = "CEMC_hotTowers_status";

  // ---- Photon builder --------------------------------------------------

  double photonMinEtGeV = 5.0;

  bool photonUseVertexCut = false;
  double photonVertexCutCm = 60.0;

  double showerShapeTowerMinEnergyGeV = 0.070;

  bool useTopoclusterIsolation = false;
  std::string topoclusterNode = "TOPOCLUSTER_ALLCALO";

  bool useSubtractedIsolation = false;

  // ---- Au+Au underlying-event subtraction --------------------------------

  double retowerFractionCut = 0.5;

  // 0 no flow, 1 psi2 from the calorimeter, 2 from HIJING, 3 from the sEPD.
  int flowModulation = 0;

  int firstPassSeedJetD = 3;
  double secondPassSeedPtGeV = 7.0;

  // ---- Archived p+p double-interaction reconstruction ------------------
  //
  // The archived Geant4-only inputs were produced with a fixed sequence of
  // random seeds consumed in order. Reproducing the donor requires replaying
  // that sequence exactly, and refusing any other seed source.

  std::vector<unsigned int> historicalSeeds;
  int expectedPedestalSequence = -1;
  bool forbidRecoConstsRandomSeed = true;
};

// The fully resolved production job. After LoadPlan and Validate the macro
// needs no knowledge of the YAML.
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

//
// The production interface, in the order the macro calls it.
//
// Parse the configuration file, resolve one named profile and the physical
// input list, and build the typed producer configuration. Calls Validate.
Plan LoadPlan(const std::string& yamlFile, const Job& job);

// Reject an internally inconsistent plan before any file is opened.
void Validate(const Plan& plan);

// Bind the run context, resolve conditions payloads through the database and
// verify every pinned digest. Updates the resolved payload paths in the plan.
void ConfigureConditions(Plan& plan);

// Register, in order: the exposure observer, the reconstruction chain the
// profile needs, and finally the producer itself.
void RegisterReconstruction(Fun4AllServer* server,
                            const Plan& plan,
                            PhotonJetTree* producer);

// Register the physical input streams and, when the calorimeter is
// reconstructed during this pass, the geometry run node.
void RegisterInputs(Fun4AllServer* server, const Plan& plan);

// Concise startup report of the resolved job.
void Print(const Plan& plan);

}  // namespace production
}  // namespace photonjet

#endif  // PHOTONJETTREE_INTERNAL_PRODUCTION_H
