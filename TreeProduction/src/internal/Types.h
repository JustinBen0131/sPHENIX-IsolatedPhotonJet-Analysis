#ifndef PHOTONJETTREE_INTERNAL_TYPES_H
#define PHOTONJETTREE_INTERNAL_TYPES_H

//
// Types.h
//
// The one typed vocabulary of the isolated-photon + jet tree producer.
//
// Everything the producer configures, measures, or writes is declared here as
// a plain value type. There is exactly one definition of each concept; the
// producer class and its implementation files use these types and never
// redeclare them.
//
// This file holds data only: enums, identities, configuration, records, and
// the pure integer functions that build relational identities. It never
// touches Fun4All, PHCompositeNode, ROOT, YAML, detector reconstruction, or a
// classification model.
//
//
// Conventions used throughout
//
//
// Units
//   energy, transverse energy, transverse momentum   GeV
//   angles                                           radians
//   positions                                        cm
//   calorimeter time                                 ADC samples, or ns where
//                                                    a field name says Ns
//   centrality                                       percent, 0 = most central
//
// Numeric types
//   physical observables    double
//   states, codes, counts   fixed-width integers
//   identities and masks    64-bit unsigned
//
//   Upstream sPHENIX accessors return float for most calorimeter quantities.
//   Widening to double is lossless and gives every observable one rule.
//
// Missing values
//   A physical observable that was not measured is NaN. A state that was not
//   determined is its Unknown enumerator, never a silently plausible default.
//   Absent, invalid, out-of-domain and not-applicable stay distinguishable,
//   because efficiency and denominator work downstream needs them separate.
//
// Identity scope
//   Identity is a compact 128-bit key joining rows of one output file. It is
//   built from the source content hash and source-local coordinates, so it is
//   unique within a source and carries no meaning across unrelated sources.
//   Full provenance (hashes, run, segment, entry ordinals) stays stored
//   explicitly and is never replaced by the compact key. A cross-file join
//   must carry the source identity with it.
//
// Ordering
//   A field named "ordinal" is encounter order: the order in which the
//   producer met the object while reading its container. A field named
//   "order" is the deterministic retained order defined by the owning capture
//   step. Neither is a transverse-momentum ranking.
//
#include <array>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace photonjet
{

//
// 1. Production vocabulary
//
// These enumerators are persisted. Their integer values are part of the output
// schema; extend the lists rather than renumbering them.
//
enum class CollisionSystem : std::uint8_t
{
  PP   = 1,
  AuAu = 2
};

enum class DataKind : std::uint8_t
{
  Data       = 1,
  Simulation = 2
};

// Which simulated sample produced this file. This names the generator lane,
// not the role a given truth object plays in the analysis signal definition.
enum class SimulationRole : std::uint8_t
{
  None         = 0,
  PhotonJet    = 1,
  InclusiveJet = 2,
  MinimumBias  = 3
};

// ReconstructedDST  the input already carries calibrated towers and clusters.
// ArchivedG4Only    the input carries Geant4 hits only, so the calorimeter and
//                   the event geometry are reconstructed during this pass.
enum class InputMode : std::uint8_t
{
  ReconstructedDST = 1,
  ArchivedG4Only   = 2
};

// Physically distinct reconstructed jet collections. In Au+Au the underlying
// event is subtracted for Sub1 and left in place for NoSub. Both are retained
// so that the subtraction itself stays measurable.
enum class JetView : std::uint8_t
{
  PP        = 1,
  AuAuSub1  = 2,
  AuAuNoSub = 3
};

// How an isolation cone sum was formed.
enum class IsolationMethod : std::int32_t
{
  Unknown         = 0,
  CalorimeterRaw  = 1,  // calibrated towers, no underlying-event subtraction
  CalorimeterSub1 = 2,  // Au+Au underlying-event subtracted towers
  Topocluster     = 3   // topological clusters, candidate removed explicitly
};

// A three-valued fact, for cases where an absent input and a measured negative
// answer must not collapse onto the same stored value.
enum class TriState : std::int32_t
{
  Unknown = -1,
  False   = 0,
  True    = 1
};

// Progress of a capture step that can legitimately do nothing.
enum class CaptureState : std::int32_t
{
  InputUnavailable = -1,  // the step ran and its required input was missing
  NotApplicable    = 0,   // the step does not apply to this production mode
  Complete         = 1    // the step ran and its input was present
};

// Relation between a reconstructed object and a truth object.
//
//   Unmatched  both exist, this pair was considered and not selected
//   Fake       reconstructed object with no corresponding truth object
//   Miss       truth object with no corresponding reconstructed object
//
// Miss may only be asserted when reconstructed capture for the event is known
// to be complete. Otherwise the relation stays Unknown, because an incomplete
// reconstruction cannot distinguish a real inefficiency from a bookkeeping gap.
enum class AssociationState : std::int32_t
{
  Unknown       = 0,
  NotApplicable = 1,
  Unmatched     = 2,
  Matched       = 3,
  Fake          = 4,
  Miss          = 5
};

// State of the calorimeter truth evaluator for one reconstructed cluster. It
// describes the evaluator, not the final reconstruction-to-truth relation.
//
// The integer values follow the reference production so that the two remain
// directly comparable during migration.
enum class DominantTruthState : std::int32_t
{
  NotSimulation        = -1,
  EvaluatorUnavailable = 0,
  NoPrimary            = 1,
  ValidPrimary         = 2,
  InvalidPrimary       = 3
};

enum class DominantTruthEvaluator : std::int32_t
{
  NotApplicable  = -1,
  Unavailable    = 0,
  TowerInfo      = 1,
  LegacyRawTower = 2
};

// Generator ancestry class of a truth photon, from the production-vertex walk.
//
//   Unknown        no generator association, or no usable production history
//   Unclassified   a production vertex was found but matched no pattern
//   Direct         two incoming and two outgoing partonic or electroweak legs
//   Fragmentation  one incoming leg radiates a photon and keeps its flavour
//   Hadronic       decay of a hadronic parent
enum class TruthPhotonClass : std::int32_t
{
  Unknown       = -1,
  Unclassified  = 0,
  Direct        = 1,
  Fragmentation = 2,
  Hadronic      = 3
};

// The truth signal predicate admits every class strictly below Hadronic. That
// deliberately includes Unknown and Unclassified: a photon whose generator
// history could not be read is not thereby a decay photon, and excluding it
// would bias the signal denominator.
inline bool isAnalysisSignalClass(const TruthPhotonClass photonClass)
{
  return static_cast<std::int32_t>(photonClass) <
         static_cast<std::int32_t>(TruthPhotonClass::Hadronic);
}

// The nominal away-side witness for a photon-jet pair, in radians. The
// comparison is inclusive: a pair is a recoil witness when its absolute
// azimuthal separation is at least this value.
inline constexpr double kRecoilDeltaPhiMin =
    7.0 * 3.14159265358979323846 / 8.0;

// Output.cc writes these enumerators as fixed-width ROOT leaves.
static_assert(sizeof(CollisionSystem) == sizeof(std::uint8_t),
              "CollisionSystem is persisted as an 8-bit value");
static_assert(sizeof(JetView) == sizeof(std::uint8_t),
              "JetView is persisted as an 8-bit value");
static_assert(sizeof(AssociationState) == sizeof(std::int32_t),
              "AssociationState is persisted as a 32-bit value");
static_assert(sizeof(DominantTruthState) == sizeof(std::int32_t),
              "DominantTruthState is persisted as a 32-bit value");
static_assert(sizeof(TriState) == sizeof(std::int32_t),
              "TriState is persisted as a 32-bit value");

//
// 2. Relational identity
//
// A compact key joining rows of one output file. (0,0) is the null identity
// and never names a real object.
struct Identity
{
  std::uint64_t hi = 0;
  std::uint64_t lo = 0;

  bool valid() const
  {
    return hi != 0 || lo != 0;
  }
};

inline bool operator==(const Identity& lhs, const Identity& rhs)
{
  return lhs.hi == rhs.hi && lhs.lo == rhs.lo;
}

inline bool operator!=(const Identity& lhs, const Identity& rhs)
{
  return !(lhs == rhs);
}

// Deterministic 64-bit mixing. It spreads already meaningful integers over the
// key space so joins do not collide. It is not a cryptographic hash and it
// carries no physics.
inline std::uint64_t mixIdentityWord(std::uint64_t value)
{
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

// Read 16 hexadecimal characters starting at offset. Returns 0 when the text
// is too short or holds a non-hexadecimal character; callers treat that as an
// unusable source identity rather than as the value zero.
inline std::uint64_t identityWordFromHex(const std::string& text,
                                         const std::size_t offset)
{
  if (text.size() < offset + 16)
  {
    return 0;
  }

  std::uint64_t value = 0;

  for (std::size_t i = 0; i < 16; ++i)
  {
    const char character = text[offset + i];

    int digit = -1;
    if (character >= '0' && character <= '9')
    {
      digit = character - '0';
    }
    else if (character >= 'a' && character <= 'f')
    {
      digit = 10 + (character - 'a');
    }
    else if (character >= 'A' && character <= 'F')
    {
      digit = 10 + (character - 'A');
    }

    if (digit < 0)
    {
      return 0;
    }

    value = (value << 4U) | static_cast<std::uint64_t>(digit);
  }

  return value;
}

// The source identity is the leading 128 bits of the pinned SHA256 of the
// input file, or of the source manifest when a single file hash does not
// apply. The full digest stays stored in SourceRecord.
inline Identity makeSourceIdentity(const std::string& sha256Hex)
{
  Identity id;
  id.hi = identityWordFromHex(sha256Hex, 0);
  id.lo = identityWordFromHex(sha256Hex, 16);
  return id;
}

// The event key binds the source to the original source entry and records
// separately whether the physical event sequence was readable, so that an
// event with an unknown sequence cannot collide with one whose sequence is
// zero.
inline Identity makeEventIdentity(const Identity& source,
                                  const std::int64_t sourceEntry,
                                  const std::int64_t physicalEventSequence,
                                  const bool physicalEventSequenceValid)
{
  const std::uint64_t entryBits =
      static_cast<std::uint64_t>(sourceEntry);

  const std::uint64_t sequenceBits =
      physicalEventSequenceValid
          ? static_cast<std::uint64_t>(physicalEventSequence)
          : 0xffffffffffffffffULL;

  const std::uint64_t sequenceTag =
      physicalEventSequenceValid ? 0x71d67fffeda60001ULL
                                 : 0x71d67fffeda60000ULL;

  Identity id;
  id.hi = mixIdentityWord(source.hi ^ entryBits);
  id.lo = mixIdentityWord(source.lo ^ sequenceBits ^ sequenceTag);
  return id;
}

inline Identity makePhotonIdentity(const Identity& event,
                                   const std::uint64_t nativeClusterKey,
                                   const std::uint32_t encounterOrdinal)
{
  Identity id;
  id.hi = mixIdentityWord(event.hi ^ nativeClusterKey);
  id.lo = mixIdentityWord(event.lo ^ (nativeClusterKey << 1U) ^
                          static_cast<std::uint64_t>(encounterOrdinal));
  return id;
}

// Radius enters a key as an integer in units of 0.01 so that the key never
// depends on a binary floating-point representation.
inline std::uint64_t radiusCode(const double radius)
{
  const double scaled = radius * 100.0;
  const double rounded = scaled < 0.0 ? scaled - 0.5 : scaled + 0.5;
  return static_cast<std::uint64_t>(static_cast<std::int64_t>(rounded));
}

inline Identity makeJetIdentity(const Identity& event,
                                const JetView view,
                                const double radius,
                                const std::uint64_t nativeKey,
                                const std::uint32_t deterministicOrder)
{
  const std::uint64_t collection =
      (static_cast<std::uint64_t>(view) << 32U) ^ radiusCode(radius);

  Identity id;
  id.hi = mixIdentityWord(event.hi ^ collection ^ nativeKey);
  id.lo = mixIdentityWord(event.lo ^ (collection << 1U) ^
                          static_cast<std::uint64_t>(deterministicOrder));
  return id;
}

inline Identity makeTruthPhotonIdentity(const Identity& event,
                                        const int embeddingId,
                                        const int nativeTrackId)
{
  Identity id;
  id.hi = mixIdentityWord(
      event.hi ^
      static_cast<std::uint64_t>(static_cast<std::int64_t>(embeddingId)));
  id.lo = mixIdentityWord(
      event.lo ^
      static_cast<std::uint64_t>(static_cast<std::int64_t>(nativeTrackId)));
  return id;
}

inline Identity makeTruthJetIdentity(const Identity& event,
                                     const double radius,
                                     const std::uint64_t nativeKey)
{
  Identity id;
  id.hi = mixIdentityWord(event.hi ^ radiusCode(radius));
  id.lo = mixIdentityWord(event.lo ^ nativeKey);
  return id;
}

inline Identity makePairIdentity(const Identity& photon,
                                 const Identity& jet,
                                 const std::uint32_t jetRank)
{
  Identity id;
  id.hi = mixIdentityWord(photon.hi ^ jet.hi);
  id.lo = mixIdentityWord(photon.lo ^ jet.lo ^
                          static_cast<std::uint64_t>(jetRank));
  return id;
}

// Relation families that could otherwise produce the same key from the same
// pair of objects.
enum class LinkKind : std::uint64_t
{
  PhotonToTruthPhoton = 1,
  JetToTruthJet       = 2
};

inline Identity makeLinkIdentity(const Identity& event,
                                 const Identity& reconstructed,
                                 const Identity& truth,
                                 const LinkKind kind)
{
  const std::uint64_t kindBits = static_cast<std::uint64_t>(kind);

  Identity id;
  id.hi = mixIdentityWord(event.hi ^ reconstructed.hi ^ truth.hi ^ kindBits);
  id.lo = mixIdentityWord(event.lo ^ reconstructed.lo ^ truth.lo ^
                          (kindBits << 32U));
  return id;
}

//
// 3. Configuration
//
// One resolved Config describes one production job. The production layer builds
// it from the shared configuration file plus per-job values, and the producer
// treats it as read-only.
//
// Translates a physical event sequence back to the original source entry.
//
// Arguments: run, physical event sequence, whether that sequence is valid.
// Returns:   the source entry ordinal of the event.
//
// Modes that read paired inputs, or that skip into the middle of a file, need
// this because the producer's own encounter count is not the original cursor.
// Modes that read one file sequentially leave it empty, and the producer uses
// firstEntry plus the encounter count.
using SourceEntryResolver =
    std::function<std::int64_t(int, std::int64_t, bool)>;

// Which physical input this job reads, and how to name it afterwards.
struct SourceConfig
{
  std::string dataset;   // production lane token
  std::string sample;    // generator sample, or data period sample name
  std::string period;    // running period this source belongs to
  std::string siDiRole;  // single- or double-interaction role, where it applies

  std::string inputUriSha256;        // digest of the resolved input location
  std::string inputFileSha256;       // digest of the input file contents
  std::string sourceManifestSha256;  // digest of the manifest listing it

  int run = 0;
  int segment = 0;

  std::int64_t sourceFileOrdinal = -1;  // position of this file in the job

  std::int64_t firstEntry = 0;  // first source entry this job reads
  std::int64_t maxEvents = -1;  // event budget, negative means no limit

  SourceEntryResolver sourceEntryResolver;
};

// Names of the input nodes the producer reads. They are configurable because
// embedded simulation and archived reconstruction place the same physical
// quantity under different node names.
struct NodeConfig
{
  std::string eventHeader = "EventHeader";
  std::string gl1Packet = "GL1Packet";

  std::string mbdVertexMap = "MbdVertexMap";
  std::string globalVertexMap = "GlobalVertexMap";
  std::string mbdGeometry = "MbdGeom";

  std::string centralityInfo = "CentralityInfo";
  std::string minimumBiasInfo = "MinimumBiasInfo";

  std::string mbdOut = "MbdOut";
  std::string mbdPmts = "MbdPmtContainer";

  std::string cemcTowers = "TOWERINFO_CALIB_CEMC";
  std::string ihcalTowers = "TOWERINFO_CALIB_HCALIN";
  std::string ohcalTowers = "TOWERINFO_CALIB_HCALOUT";

  std::string cemcGeometry = "TOWERGEOM_CEMC";

  std::string truthInfo = "G4TruthInfo";
  std::string hepmcEventMap = "PHHepMCGenEventMap";
};

struct CalorimeterConfig
{
  // Whether event energy sums admit only towers the calibration marked good.
  // The choice is recorded in every event so that a reader never has to guess
  // which convention produced a stored sum.
  bool requireGoodTowers = true;
  int backgroundFlowMode = 0;  // production steering, recorded for AuAu only

};

// Photon capture.
//
// These bounds are the storage domain, not the measurement selection. The
// analysis window is narrower and lives downstream. Storing the wider domain
// is what lets a working point move without returning to the DST.
struct PhotonConfig
{
  double minEtGeV = 5.0;
  double maxEtGeV = 40.0;
  double maxAbsEta = 0.7;

  // Reconstructed-object vertex domain, cm. Photons and jets are captured
  // only when the reconstructed vertex is valid and |z| is below this value;
  // the event row is retained either way and records the outcome. This is
  // the explicit object domain of the retention profile, not an analysis
  // vertex selection, which is narrower and lives downstream.
  double objectVertexAbsZMaxCm = 60.0;

  std::string rawClusterNode = "CLUSTERINFO_CEMC";
  std::string photonNode = "PHOTONCLUSTER_CEMC";
  std::string topoclusterNode = "TOPOCLUSTER_ALLCALO";

  // Reconstructed isolation cone radii retained for every candidate.
  std::vector<double> isolationRadii = {0.30, 0.40};

  // Truth isolation. The signal predicate uses the signal radius.
  double truthIsolationSignalRadius = 0.30;
  double truthIsolationMaxEtGeV = 4.0;
  double truthIsolationCoreRadius = 0.001;

  // Conversion from calorimeter time samples to nanoseconds.
  double timingSampleNs = 17.6;
};

// One reconstructed jet collection: a radius, a physical view, and the pair of
// nodes holding the same jets before and after the energy-scale correction.
struct JetNodeConfig
{
  JetView view = JetView::PP;
  double radius = 0.4;

  std::string rawNode;        // before the jet energy-scale correction
  std::string correctedNode;  // after it
  std::string truthNode;      // matching truth collection, simulation only

  std::string inputIdentity;        // tower input description
  std::string subtractionIdentity;  // underlying-event treatment description
};

struct JetConfig
{
  std::vector<JetNodeConfig> nodes;

  // Storage domain, not the recoil selection.

  // Geometric radius within which a reconstructed jet may be associated with a
  // truth jet of the same clustering radius.
  double truthMatchMaxDeltaR = 0.30;

  // The jet energy scale is applied once, upstream, by the reconstruction the
  // production layer registers. The producer stores the payload identity so a
  // reader can tell which correction produced the stored corrected momentum.
  std::string energyScalePayloadPath;
  std::string energyScalePayloadSha256;
  bool energyScaleUsesEmFraction = false;
};

struct CentralityConfig
{
  // Coarse analysis bin edges in percent, stored beside the native centile as
  // a join convenience. The native value stays authoritative.
  std::vector<int> edges = {0, 20, 50, 80};

  // Au+Au data reconstructs centrality during this pass from run-bound frozen
  // payloads. Embedded simulation consumes the centrality already present in
  // its input instead.
  bool reconstructDataCentrality = false;
  bool requireOwnRunPayload = true;

  std::string calibrationTag;
  std::string divisionPayloadPath;
  std::string runScalePayloadPath;
  std::string vertexScalePayloadPath;
};

// Producer-side weight inputs.
//
// The only factor the producer itself applies is the simulated-vertex
// reweighting, a lookup in a histogram of data-over-simulation vertex
// probability keyed by the truth hard-scatter vertex. Slice cross sections,
// sample stitching, interaction mixing and luminosity are downstream, where
// they are combined once from the sample manifest.
struct WeightConfig
{
  bool applyVertexReweight = false;

  std::string vertexReweightFile;
  std::string vertexReweightFileSha256;
  std::string vertexReweightHistogram;
};

// Which optional tables the file carries. The mandatory tables are always
// written; these switches govern the large witness tables whose final
// disposition is a mapping-ledger decision.
struct OutputConfig
{
  bool writePhotonCells = true;
  bool writeIsolationConstituents = true;
  bool writePhotonJetPairs = true;
};

// Everything needed to say which code, configuration and calibration produced
// a file. Written once into the output metadata.
struct ProvenanceConfig
{
  std::string productionTag;
  std::string softwareRelease;

  std::string producerGitCommit;
  std::string producerSourceSha256;
  std::string producerLibrarySha256;

  std::string macroSha256;
  std::string configurationSha256;

  // The exact configuration bytes, embedded in the output so a file never
  // has to be matched to a configuration by name.
  std::string configurationText;

  std::string calibrationManifestSha256;
};

struct Config
{
  CollisionSystem system = CollisionSystem::PP;
  DataKind dataKind = DataKind::Data;
  SimulationRole simulationRole = SimulationRole::None;
  InputMode inputMode = InputMode::ReconstructedDST;

  std::string outputFile;

  SourceConfig source;
  NodeConfig nodes;
  CalorimeterConfig calorimeter;
  PhotonConfig photon;
  JetConfig jets;
  CentralityConfig centrality;
  WeightConfig weights;
  OutputConfig output;
  ProvenanceConfig provenance;

  // Names of the shower-shape definitions to evaluate for every candidate.
  // "H70" and "H0" are the hybrid definition (rectangular sums over the full
  // good-tower grid, moments over cluster-owned towers) at a 70 MeV and a 0 MeV
  // floor; these are the two definitions the identification models have used.
  std::vector<std::string> showerDefinitions = {"H70", "H0"};

  // Write an event row for every event the producer is handed, including
  // events with no photon candidate. The event population is the exposure
  // denominator of every later normalisation, so it is not filtered by object
  // yield.
  bool retainEveryProducerEncounter = true;

  // Embedded Au+Au simulation may require the embedded minimum-bias
  // classifier. The event record is retained either way; the requirement
  // affects object capture, not exposure accounting.
  bool requireEmbeddedMinimumBias = false;
};

//
// 4. Source and event accounting
//
// One row per input source, written once at the end when the counts are final.
// A reader uses it to confirm that a file covers the entry range it claims and
// that the job finished rather than stopped.
struct SourceRecord
{
  Identity sourceId;

  std::string dataset;
  std::string sample;
  std::string period;
  std::string siDiRole;

  std::string inputUriSha256;
  std::string inputFileSha256;
  std::string sourceManifestSha256;

  int run = 0;
  int segment = 0;

  std::int64_t sourceFileOrdinal = -1;

  std::int64_t firstEntry = 0;
  std::int64_t lastEntry = -1;

  // encounteredEvents       events this producer was handed
  // retainedEvents          events written to the event tree
  // upstreamRejectedEvents  events present in the input but stopped by
  //                         reconstruction before reaching this producer
  //
  // The exposure denominator is encountered plus upstream rejected. The
  // producer cannot observe rejected events by itself; see
  // UpstreamRejectedEventRecord.
  std::uint64_t encounteredEvents = 0;
  std::uint64_t retainedEvents = 0;
  std::uint64_t upstreamRejectedEvents = 0;

  bool completed = false;
};

// One row per input event that reconstruction stopped before this producer ran.
//
// A module in the reconstruction chain can abort an event, and Fun4All then
// skips every module after it, including this one. Those events are real
// exposure and must not vanish from the denominator. They are supplied by a
// companion observer registered ahead of reconstruction, through
// PhotonJetTree::recordUpstreamRejectedEvent.
struct UpstreamRejectedEventRecord
{
  Identity sourceId;

  int run = 0;
  std::int64_t physicalEventSequence = -1;
  std::int64_t sourceEntry = -1;
};

// One row per retained event.
struct EventRecord
{
  Identity sourceId;
  Identity eventId;

  //
  // Where this event came from
  //
  // The physical event sequence is what the data acquisition assigned. The
  // producer event ordinal is this job's own counter. The source entry is the
  // position in the original file. All three are stored because a job that
  // skips into a file has all three differ.
  //
  int run = 0;
  int segment = 0;

  std::int64_t sourceFileOrdinal = -1;
  std::int64_t sourceEntry = -1;
  std::int64_t sourceGlobalEntry = -1;

  std::uint64_t producerEventOrdinal = 0;

  std::int64_t physicalEventSequence = -1;
  bool physicalEventSequenceValid = false;

  //
  // Trigger
  //
  // The three decision words are read through their own accessors and are not
  // interchangeable. The availability mask records which of them the packet
  // version actually supplies, so a zero that was never measured cannot be
  // mistaken for a measured zero.
  //
  // Historical getTriggerVector witness: aliases the live word in v2/v3.
  // It must not be interpreted as an independently measured raw decision.
  std::uint64_t triggerInputBits = 0;
  std::uint64_t triggerLiveBits = 0;
  std::uint64_t triggerScaledBits = 0;

  std::uint32_t triggerPacketNumber = 0;
  std::uint64_t triggerBunchNumber = 0;
  std::uint64_t triggerPacketStatus = 0;

  std::int32_t triggerPacketVersion = 0;
  std::uint32_t triggerDecisionsAvailable = 0;
  bool triggerPacketValid = false;

  // Index into the scaler snapshot table. Consecutive events that saw no
  // counter change share one snapshot.
  std::int64_t triggerScalerSnapshotId = -1;

  //
  // Reconstructed vertex
  //
  // Validity is independent of the stored number: an invalid vertex stores NaN
  // rather than a plausible zero.
  //
  double recoVertexX = std::numeric_limits<double>::quiet_NaN();
  double recoVertexY = std::numeric_limits<double>::quiet_NaN();
  double recoVertexZ = std::numeric_limits<double>::quiet_NaN();
  bool recoVertexValid = false;

  // 0 none, 1 the minimum-bias detector vertex map, 2 the global vertex map.
  std::int32_t recoVertexSource = 0;

  // True when the reconstructed vertex is valid and inside the configured
  // object vertex domain, so that photon and jet capture ran for this event.
  bool recoObjectVertexInDomain = false;

  //
  // Minimum-bias detector
  //
  // The per-tube arrays are the inputs of the centrality calibration. They are
  // kept so that centrality can be recomputed from a later calibration without
  // returning to the DST. The charge population the centrality calibration
  // selects is not the same thing as a positive-charge sum.
  //
  double mbdT0Ns = std::numeric_limits<double>::quiet_NaN();
  double mbdSouthTimeNs = std::numeric_limits<double>::quiet_NaN();
  double mbdNorthTimeNs = std::numeric_limits<double>::quiet_NaN();

  double mbdSouthCharge = std::numeric_limits<double>::quiet_NaN();
  double mbdNorthCharge = std::numeric_limits<double>::quiet_NaN();
  double mbdTotalCharge = std::numeric_limits<double>::quiet_NaN();

  bool mbdValid = false;

  TriState mbdPmtAvailable = TriState::Unknown;

  std::vector<std::int32_t> mbdPmtId;
  std::vector<std::int32_t> mbdPmtArm;
  std::vector<double> mbdPmtCharge;
  std::vector<double> mbdPmtTimeNs;
  std::vector<std::int32_t> mbdPmtValid;

  //
  // Centrality
  //
  // The value written by the reconstruction is kept separately from the
  // canonical value this producer reports, so a later recalibration never
  // destroys the original measurement.
  //
  bool centralityApplicable = false;

  TriState minimumBiasDecision = TriState::Unknown;

  double centralitySelectedCharge = std::numeric_limits<double>::quiet_NaN();

  std::int32_t centralityNativeBin = -1;
  double centralityNativePercent = std::numeric_limits<double>::quiet_NaN();
  bool centralityNativeValid = false;

  std::int32_t centralityBin = -1;
  double centralityPercent = std::numeric_limits<double>::quiet_NaN();
  bool centralityValid = false;

  //
  // Calorimeter event state
  //
  // Sums are signed and admit finite entries only. A container that is present
  // but partly unreadable keeps its finite partial sum and loses its valid
  // bit, because the partial sum is still informative and the reader must know
  // that it is one.
  //
  double cemcEnergy = std::numeric_limits<double>::quiet_NaN();
  double ihcalEnergy = std::numeric_limits<double>::quiet_NaN();
  double ohcalEnergy = std::numeric_limits<double>::quiet_NaN();
  double totalCaloEnergy = std::numeric_limits<double>::quiet_NaN();

  // bits 0,1,2   CEMC, inner and outer hadronic container present
  // bits 8,9,10  the same containers complete and finite
  std::uint32_t caloAvailableMask = 0;
  std::uint32_t caloValidMask = 0;

  bool caloRequiredGoodTowers = true;

  // Per-detector tower census, indexed CEMC, inner hadronic, outer hadronic.
  std::array<std::uint64_t, 3> caloTowerCount{};
  std::array<std::uint64_t, 3> caloAcceptedTowerCount{};
  std::array<std::uint64_t, 3> caloBadQualityTowerCount{};
  std::array<std::uint64_t, 3> caloNullTowerCount{};
  std::array<std::uint64_t, 3> caloNonFiniteTowerCount{};

  //
  // AuAu background primitives, from TowerInfoBackground_Sub2. The three
  // vectors are detector eta strips; retain them to reproduce subtraction.
  CaptureState jetBackgroundState = CaptureState::NotApplicable;
  std::int32_t jetBackgroundFlowMode = -1;
  double jetBackgroundV2 = std::numeric_limits<double>::quiet_NaN();
  double jetBackgroundPsi2 = std::numeric_limits<double>::quiet_NaN();
  std::vector<float> jetBackgroundUEEmcal, jetBackgroundUEIhcal, jetBackgroundUEOhcal;
  std::int32_t jetBackgroundNStrips = -1, jetBackgroundNTowers = -1;
  TriState jetBackgroundFlowFailure = TriState::Unknown;

  // One entry per configured view/radius, including empty and skipped views.
  std::vector<std::int32_t> recoJetView, recoJetRadiusCode, recoJetCaptureState;
  std::vector<std::uint32_t> recoJetCount;

  // Simulation truth context
  //
  CaptureState truthVertexCaptureState = CaptureState::NotApplicable;

  std::int32_t truthPrimaryVertexId = -1;

  double truthHardVertexZ = std::numeric_limits<double>::quiet_NaN();
  bool truthHardVertexValid = false;

  double truthMinimumBiasVertexZ = std::numeric_limits<double>::quiet_NaN();
  bool truthMinimumBiasVertexValid = false;

  TriState embeddedMinimumBias = TriState::Unknown;

  // Truth photon census. These counters are what makes the truth denominator
  // auditable: a reader can tell how many candidate truth photons existed, how
  // many were written, and why the rest were not.
  CaptureState truthPhotonCaptureState = CaptureState::NotApplicable;

  std::uint32_t truthPhotonEmbeddedPrimaryCount = 0;
  std::uint32_t truthPhotonWrittenCount = 0;
  std::uint32_t truthPhotonRejectedCount = 0;
  std::uint32_t truthPhotonRejectedKinematicsCount = 0;
  std::uint32_t truthPhotonRejectedIsolationCount = 0;
  std::uint32_t truthPhotonDuplicateTrackCount = 0;
  std::uint32_t truthPhotonIsolationIncompleteCount = 0;
  std::uint32_t truthPhotonAnalysisSignalCount = 0;

  // False when duplicate identities or incomplete isolation inputs make the
  // truth-side miss bookkeeping unreliable for this event.
  bool truthDenominatorComplete = false;

  CaptureState truthJetCaptureState = CaptureState::NotApplicable;

  std::vector<std::int32_t> truthJetRadiusCode;
  std::vector<std::int32_t> truthJetContainerValid;

  //
  // Reconstructed photon completeness
  //
  // A truth photon may only be called a miss when reconstructed capture is
  // known to be complete for the event.
  //
  CaptureState recoPhotonCaptureState = CaptureState::NotApplicable;
  std::uint32_t recoPhotonUnclassifiableCount = 0;

  //
  // Producer weight
  //
  // One composed weight per event, with its factors in the weight table. Data
  // carries exactly one, which is a producer convention and not a luminosity
  // normalisation. Cross-sample normalisation is downstream.
  //
  double eventWeight = std::numeric_limits<double>::quiet_NaN();
  bool eventWeightValid = false;

  //
  // Object census for this event
  //
  std::uint32_t photonCount = 0;
  std::uint32_t jetCount = 0;
  std::uint32_t photonJetPairCount = 0;
  std::uint32_t truthPhotonCount = 0;
  std::uint32_t truthJetCount = 0;

  // Producer terminal status. Zero is ordinary completion. Non-zero values are
  // diagnostic and never act as an output retention cut.
  std::int32_t terminalStatus = 0;
};

//
// 5. Reconstructed photons
//
// Cluster timing.
//
// The mean is energy weighted over the towers the cluster owns, with no
// quality gate and no energy floor, which is the definition the reconstruction
// uses. Numerator, denominator and contributing count are stored so that the
// value is reproducible and so that a sentinel can never be mistaken for a
// measurement.
struct PhotonTimingRecord
{
  double meanTimeSamples = std::numeric_limits<double>::quiet_NaN();

  double energyWeightedNumerator = std::numeric_limits<double>::quiet_NaN();
  double energyDenominator = std::numeric_limits<double>::quiet_NaN();

  std::uint32_t contributingTowerCount = 0;

  bool finite = false;
  bool valid = false;

  // meanTimeSamples converted with PhotonConfig::timingSampleNs.
  double timeNs = std::numeric_limits<double>::quiet_NaN();
};

// One shower-shape definition evaluated for one photon candidate.
//
// The complete rectangular family is stored, not a chosen subset, because the
// identification models are ratios and widths over these sums. With the family
// present, together with the event vertex and the centrality, the registered models' proposed
// input vector is reconstructible downstream without returning to the DST.
// That is what makes the base tree model-independent in substance rather than
// only in name.
//
// Naming follows the reconstruction: eNM is the rectangular sum over N towers
// in pseudorapidity by M in azimuth about the seed tower. e32 is the
// three-by-two block on the side of the azimuthal centre of gravity.
//
// Widths are energy-weighted second moments about the centre of gravity, with
// the central tower excluded from the numerator and included in the
// denominator. No square root is taken.
//
// A candidate may carry more than one definition when a production needs to
// compare them; definitionName says which one this row is.
struct ShowerShapeRecord
{
  std::string definitionName;

  // Tower energy floor of this definition, GeV. Towers at or below it do not
  // enter the sums.
  double energyFloorGeV = 0.0;

  // Rectangular sums, GeV.
  double e11 = std::numeric_limits<double>::quiet_NaN();
  double e13 = std::numeric_limits<double>::quiet_NaN();
  double e15 = std::numeric_limits<double>::quiet_NaN();
  double e17 = std::numeric_limits<double>::quiet_NaN();
  double e22 = std::numeric_limits<double>::quiet_NaN();
  double e31 = std::numeric_limits<double>::quiet_NaN();
  double e32 = std::numeric_limits<double>::quiet_NaN();
  double e33 = std::numeric_limits<double>::quiet_NaN();
  double e35 = std::numeric_limits<double>::quiet_NaN();
  double e37 = std::numeric_limits<double>::quiet_NaN();
  double e51 = std::numeric_limits<double>::quiet_NaN();
  double e52 = std::numeric_limits<double>::quiet_NaN();
  double e53 = std::numeric_limits<double>::quiet_NaN();
  double e55 = std::numeric_limits<double>::quiet_NaN();
  double e57 = std::numeric_limits<double>::quiet_NaN();
  double e71 = std::numeric_limits<double>::quiet_NaN();
  double e72 = std::numeric_limits<double>::quiet_NaN();
  double e73 = std::numeric_limits<double>::quiet_NaN();
  double e75 = std::numeric_limits<double>::quiet_NaN();
  double e77 = std::numeric_limits<double>::quiet_NaN();

  // Second moments over the full support.
  double weta = std::numeric_limits<double>::quiet_NaN();
  double wphi = std::numeric_limits<double>::quiet_NaN();

  // Second moments about the centre of gravity. The cogX form uses the
  // sub-tower centre of gravity and is the one the identification models take.
  double wetaCog = std::numeric_limits<double>::quiet_NaN();
  double wphiCog = std::numeric_limits<double>::quiet_NaN();
  double wetaCogX = std::numeric_limits<double>::quiet_NaN();
  double wphiCogX = std::numeric_limits<double>::quiet_NaN();

  // Second moments restricted to the three-by-three core. The canonical Au+Au
  // identification model takes these two. They remain NaN until the upstream
  // photon builder supplies them; see the upstream reconciliation step.
  double weta33CogX = std::numeric_limits<double>::quiet_NaN();
  double wphi33CogX = std::numeric_limits<double>::quiet_NaN();

  // Azimuthal widths of the three-, five- and seven-tower strips.
  double w32 = std::numeric_limits<double>::quiet_NaN();
  double w52 = std::numeric_limits<double>::quiet_NaN();
  double w72 = std::numeric_limits<double>::quiet_NaN();

  // Shower fractions from the reconstruction's own accessor. These are
  // dimensionless combinations of an oriented two-by-two block over the total
  // cluster energy, not rectangular sums with similar names.
  double et1 = std::numeric_limits<double>::quiet_NaN();
  double et2 = std::numeric_limits<double>::quiet_NaN();
  double et3 = std::numeric_limits<double>::quiet_NaN();
  double et4 = std::numeric_limits<double>::quiet_NaN();

  // Seed tower and centre of gravity in tower index space.
  std::int32_t centerEtaIndex = -1;
  std::int32_t centerPhiIndex = -1;
  double centerOfGravityEta = std::numeric_limits<double>::quiet_NaN();
  double centerOfGravityPhi = std::numeric_limits<double>::quiet_NaN();

  // Extent of the owned towers about the seed in tower units, and the radial
  // spread about the centre of gravity.
  std::int32_t deltaEtaMax = -1;
  std::int32_t deltaPhiMax = -1;
  double deltaEtaCog = std::numeric_limits<double>::quiet_NaN();
  double deltaPhiCog = std::numeric_limits<double>::quiet_NaN();
  double radialSpread = std::numeric_limits<double>::quiet_NaN();

  std::uint32_t saturatedTowerCount = 0;

  // Cell census of this definition's support.
  std::uint32_t ownedCellCount = 0;
  std::uint32_t goodCellCount = 0;
  std::uint32_t zeroCellCount = 0;
  std::uint32_t negativeCellCount = 0;
  std::uint32_t nonFiniteCellCount = 0;

  bool valid = false;
};

// One calorimeter cell in the neighbourhood of a photon candidate.
//
// Both the calibrated tower energy and the value the cluster itself recorded
// are kept, with the full tower status word rather than a reduced good flag,
// so that a later quality convention can be applied without re-running
// reconstruction.
//
// DISPOSITION PENDING (reference-production field mapping): whether the cell
// table stays in the delivered tree or becomes an optional diagnostic product
// is a mapping-ledger decision, not one taken here.
struct PhotonCellRecord
{
  Identity eventId;
  Identity photonId;

  std::uint64_t towerKey = 0;

  // 0 CEMC, 1 inner hadronic, 2 outer hadronic.
  std::int32_t subsystem = -1;

  std::int32_t etaIndex = -1;
  std::int32_t phiIndex = -1;

  double calibratedEnergy = std::numeric_limits<double>::quiet_NaN();
  double clusterMapEnergy = std::numeric_limits<double>::quiet_NaN();
  double towerTimeSamples = std::numeric_limits<double>::quiet_NaN();

  std::uint32_t towerStatus = 0;

  bool ownedByCluster = false;
  bool inShowerSupport = false;

  bool good = false;
  bool zeroEnergy = false;
  bool negativeEnergy = false;
  bool finite = false;
};

// One isolation cone for one photon candidate.
//
// The per-layer components are stored beside the total because the layers are
// combined differently in different systems, and because a later change of
// combination must not require re-reading the DST.
struct IsolationRecord
{
  Identity eventId;
  Identity photonId;

  double axisEta = std::numeric_limits<double>::quiet_NaN();
  double axisPhi = std::numeric_limits<double>::quiet_NaN();
  double radius = std::numeric_limits<double>::quiet_NaN();
  IsolationMethod method = IsolationMethod::Unknown;

  double coneSum = std::numeric_limits<double>::quiet_NaN();

  double electromagneticComponent = std::numeric_limits<double>::quiet_NaN();
  double innerHadronicComponent = std::numeric_limits<double>::quiet_NaN();
  double outerHadronicComponent = std::numeric_limits<double>::quiet_NaN();

  // Whether the candidate's own transverse energy was removed from the sum.
  bool candidateRemoved = false;

  bool valid = false;
};

// One object entering an isolation cone.
//
// DISPOSITION PENDING (reference-production field mapping): constituent-level
// storage is large and its final disposition is a mapping-ledger decision. It
// is retained now because dropping it would be irreversible and because the
// stored cone sums cannot otherwise be re-derived.
struct IsolationConstituentRecord
{
  Identity eventId;
  Identity photonId;

  double radius = std::numeric_limits<double>::quiet_NaN();

  // What kind of object this constituent is, for example a topological cluster
  // or a calorimeter tower of a named layer.
  std::string source;

  std::uint64_t nativeKey = 0;
  std::int32_t subsystem = -1;

  double deltaEta = std::numeric_limits<double>::quiet_NaN();
  double deltaPhi = std::numeric_limits<double>::quiet_NaN();
  double deltaR = std::numeric_limits<double>::quiet_NaN();

  double rawEnergy = std::numeric_limits<double>::quiet_NaN();
  double transverseEnergy = std::numeric_limits<double>::quiet_NaN();
  double subtractedEnergy = std::numeric_limits<double>::quiet_NaN();
  double subtractedTransverseEnergy = std::numeric_limits<double>::quiet_NaN();

  std::int32_t qualityState = 0;
  bool masked = false;
  std::uint32_t maskState = 0;  // bit0 raw bad, bit1 SUB1 bad
  bool candidateRemoved = false;
};

// The Geant4 primary that deposited the most calorimeter energy in a cluster.
//
// This is an identity witness, not the analysis association: it says which
// particle made the cluster, and the truth stage decides what that means.
struct DominantTruthRecord
{
  DominantTruthState state = DominantTruthState::NotSimulation;
  DominantTruthEvaluator evaluator = DominantTruthEvaluator::NotApplicable;

  std::int32_t trackId = -1;
  std::int32_t pid = 0;
  std::int32_t barcode = std::numeric_limits<std::int32_t>::min();
  std::int32_t embeddingId = 0;
  std::int32_t vertexId = -1;

  double energyContribution = std::numeric_limits<double>::quiet_NaN();
};

// One reconstructed photon candidate.
//
// Candidates are admitted on kinematics alone. No identification score and no
// isolation working point takes part in retention, which is what lets a model
// change without re-reading the DST.
struct PhotonRecord
{
  Identity eventId;
  Identity photonId;

  std::uint64_t nativeClusterKey = 0;
  std::uint32_t encounterOrdinal = 0;

  double energy = std::numeric_limits<double>::quiet_NaN();
  double et = std::numeric_limits<double>::quiet_NaN();
  double eta = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();

  bool kinematicsFinite = false;

  // The vertex the reconstruction used to turn cluster position into
  // direction. Stored because it is part of the measurement.
  double producerVertexZ = std::numeric_limits<double>::quiet_NaN();

  PhotonTimingRecord timing;
  DominantTruthRecord dominantTruth;

  std::vector<ShowerShapeRecord> showerShapes;
  std::vector<IsolationRecord> isolation;
};

//
// 6. Reconstructed jets
//
// One reconstructed jet.
//
// Raw and corrected transverse momentum are both kept, with the collection
// identity, so a reader can tell which energy-scale correction produced the
// stored value. Area is the measured jet area from the clustering, never a
// geometric substitute.
struct JetRecord
{
  Identity eventId;
  Identity jetId;

  JetView view = JetView::PP;
  double radius = 0.0;

  std::string inputIdentity;
  std::string subtractionIdentity;

  std::uint64_t nativeKey = 0;
  std::uint64_t nativeRawKey = 0;

  // Position in the container as read, before deterministic ordering.
  std::uint32_t encounterOrdinal = 0;

  // Deterministic retained order within this view and radius, starting at
  // zero. It makes the output reproducible; it is not a momentum ranking.
  std::uint32_t deterministicOrder = 0;

  double rawPt = std::numeric_limits<double>::quiet_NaN();
  double correctedPt = std::numeric_limits<double>::quiet_NaN();
  double eta = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();
  double area = std::numeric_limits<double>::quiet_NaN();

  std::uint64_t qualityBitmask = 0;

  bool calibrationValid = false;
};

// One photon and one retained jet of one collection.
//
// The recoil state is the nominal away-side witness, an absolute azimuthal
// separation of at least kRecoilDeltaPhiMin. It is a stored flag and not a
// retention cut: every pair is written, and any away-side requirement is
// reproducible downstream from the stored separation.
//
// DISPOSITION PENDING (reference-production field mapping): these rows are a
// Cartesian product and form the largest table in the file. They may become a
// downstream derivation once the ledger shows that the retained photon and jet
// rows carry enough ordering, population and validity information to rebuild
// them exactly.
struct PhotonJetPairRecord
{
  Identity eventId;
  Identity pairId;

  Identity photonId;
  Identity jetId;

  JetView jetView = JetView::PP;
  double jetRadius = 0.0;

  double deltaPhi = std::numeric_limits<double>::quiet_NaN();

  // Transverse momentum balance, corrected jet momentum over photon transverse
  // energy. NaN when the denominator is not positive.
  double xJGamma = std::numeric_limits<double>::quiet_NaN();

  bool recoilWitness = false;

  // Historical photon_rank is always zero; jet_rank is retained encounter
  // order and restarts in each view/radius. Neither is a pT ranking.
  std::uint32_t photonRank = 0;
  std::uint32_t jetRank = 0;
};

//
// 7. Simulation truth
//
struct TruthVertexRecord
{
  Identity eventId;

  std::int32_t vertexId = -1;
  std::int32_t embeddingId = 0;
  bool embeddingValid = false;

  double z = std::numeric_limits<double>::quiet_NaN();
  bool valid = false;
};

// One generated photon of the embedded primary population.
//
// Generator identity and Geant4 identity are kept separately, and a missing
// generator association is recorded as missing rather than replaced by an
// invented ancestry. Truth isolation sums the transverse energy of all
// embedded primary particles in the cone, not only photons, with the merged
// core removed.
struct TruthPhotonRecord
{
  Identity eventId;
  Identity truthPhotonId;

  std::int32_t trackId = -1;
  std::int32_t vertexId = -1;
  std::int32_t barcode = std::numeric_limits<std::int32_t>::min();
  std::int32_t embeddingId = 0;
  std::int32_t pid = 22;

  double pt = std::numeric_limits<double>::quiet_NaN();
  double eta = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();

  bool geantValid = false;
  bool generatorAssociationValid = false;

  TruthPhotonClass promptClass = TruthPhotonClass::Unknown;

  double isolationR03 = std::numeric_limits<double>::quiet_NaN();
  double isolationR04 = std::numeric_limits<double>::quiet_NaN();
  bool isolationValid = false;

  // The analysis signal predicate: a Geant4-valid photon inside the
  // pseudorapidity acceptance, with a prompt class below Hadronic, and with a
  // valid signal-radius isolation below the configured limit.
  bool analysisSignal = false;

  // Which generator lane produced this photon, and whether it is the lane the
  // analysis takes its signal from in this production mode.
  SimulationRole sampleRole = SimulationRole::None;
  bool signalSourceRole = false;
};

struct TruthJetRecord
{
  Identity eventId;
  Identity truthJetId;

  double radius = 0.0;
  std::uint64_t nativeKey = 0;

  double pt = std::numeric_limits<double>::quiet_NaN();
  double eta = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();

  bool containerValid = false;
};

//
// 8. Reconstruction to truth relations
//
// A reconstructed photon and a truth photon.
//
// The association is by particle identity: the truth particle that deposited
// the most energy in the cluster. It is not a nearest-distance match, and the
// stored separation is a diagnostic of the association rather than its
// criterion. Several clusters may legitimately point at one truth photon.
struct PhotonTruthLinkRecord
{
  Identity linkId;
  Identity eventId;

  Identity photonId;
  Identity truthPhotonId;

  AssociationState state = AssociationState::Unknown;

  double energyContribution = std::numeric_limits<double>::quiet_NaN();
  double deltaR = std::numeric_limits<double>::quiet_NaN();
};

// A reconstructed jet and a truth jet of the same clustering radius.
//
// Every candidate edge within the matching radius is retained; selectedMatch
// marks the one chosen by the deterministic one-to-one assignment. Matching
// runs independently in each reconstructed view, so one truth jet may be
// matched once per view.
struct JetTruthLinkRecord
{
  Identity linkId;
  Identity eventId;

  Identity jetId;
  Identity truthJetId;

  JetView jetView = JetView::PP;
  double jetRadius = 0.0;

  AssociationState state = AssociationState::Unknown;

  double deltaR = std::numeric_limits<double>::quiet_NaN();
  bool selectedMatch = false;
};

//
// 9. Event weights
//
// One weight factor, or the composed producer weight, for one event.
//
// A factor that was not determined is NaN with valid false. It is never
// assumed to be one: a missing factor and a factor of unity have different
// consequences downstream. applicationCount guards against a factor being
// folded in twice.
//
// Cross-sample normalisation, sample stitching and luminosity are downstream.
struct WeightComponentRecord
{
  Identity eventId;

  // "final" names the composed producer weight. Any other value names a
  // factor that contributed to it.
  std::string componentType;

  double sliceWeight = std::numeric_limits<double>::quiet_NaN();
  double crossSectionWeight = std::numeric_limits<double>::quiet_NaN();
  double vertexWeight = std::numeric_limits<double>::quiet_NaN();
  double siDiWeight = std::numeric_limits<double>::quiet_NaN();
  double periodWeight = std::numeric_limits<double>::quiet_NaN();
  double exposureWeight = std::numeric_limits<double>::quiet_NaN();

  double finalWeight = std::numeric_limits<double>::quiet_NaN();

  bool valid = false;
  std::int32_t applicationCount = 0;
};

//
// 10. Trigger bookkeeping
//
// A snapshot of the cumulative trigger scalers.
//
// The counters are cumulative and are not a luminosity. A new row is written
// only when a counter changes, and the event and bunch interval it covers is
// recorded, so consecutive unchanged events share one row. A counter that goes
// backwards sets the discontinuity flag rather than being silently accepted.
struct TriggerScalerSnapshotRecord
{
  Identity sourceId;

  std::int64_t snapshotId = -1;

  std::uint64_t bunchStart = 0;
  std::uint64_t bunchEnd = 0;

  std::int64_t sourceEntryStart = -1;
  std::int64_t sourceEntryEnd = -1;

  std::uint64_t observedEvents = 0;

  std::array<std::uint64_t, 64> raw{};
  std::array<std::uint64_t, 64> live{};
  std::array<std::uint64_t, 64> scaled{};

  bool valid = false;
  bool discontinuity = false;
};

// Run-level trigger configuration, one row per trigger bit.
struct TriggerRunInfoRecord
{
  Identity sourceId;

  int run = 0;
  std::int32_t bit = -1;

  std::string name;

  // The prescale the hardware was configured with, and the prescale the run
  // actually averaged. They differ when the configuration changed mid-run.
  double hardwarePrescale = std::numeric_limits<double>::quiet_NaN();
  double runAveragePrescale = std::numeric_limits<double>::quiet_NaN();

  std::uint64_t rawCount = 0;
  std::uint64_t liveCount = 0;
  std::uint64_t scaledCount = 0;

  bool valid = false;
};

}  // namespace photonjet

#endif  // PHOTONJETTREE_INTERNAL_TYPES_H
