#ifndef PHOTONJETTREE_INTERNAL_TYPES_H
#define PHOTONJETTREE_INTERNAL_TYPES_H

/**
 * @file Types.h
 * Shared configuration, identities and persisted record definitions.
 * Units, invalid states and indexing domains accompany the relevant fields.
 * These records are the common capture/persistence contract for all presets.
 * Energy and momentum: GeV; angles: rad; positions: cm; centrality: %.
 * Time is ADC samples unless the field explicitly says Ns. Unmeasured
 * observables are NaN; state/validity fields distinguish missing cases.
 * Ordinals mean encounter order, never momentum rank. Compact identities
 * join within a source; cross-file joins also carry full source provenance.
 */

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
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
// The integer values are the stored eventTree codes.
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
  std::string lane;      // production profile (lane) name, e.g. pp_photon_sim
  std::string dataset;   // production lane token
  std::string sample;    // generator sample, or data period sample name
  std::string period;    // running period this source belongs to
  std::string siDiRole;  // single- or double-interaction role, where it applies

  std::string inputUriSha256;        // digest of the resolved input location
  std::string inputFileSha256;       // digest of the input file contents
  std::string sourceManifestSha256;  // digest of the manifest listing it

  int run = 0;
  int segment = 0;

  // Position of this file in the maintained source order, and the global
  // entry of its entry 0 in that order (-1 when the catalogue has none). The
  // PPG12 half-sample convention keys on this global entry.
  std::int64_t sourceFileOrdinal = -1;
  std::int64_t sourceGlobalEntryBegin = -1;

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
  // The embedding fallback is only for collaborator charge capture; the
  // centrality replay remains bound to the configured primary PMT node.
  std::string mbdPmtsFallback = "MbdPmtContainer_data";

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

// One retained shower-feature view definition.
//
// energySource: 0 calibrated TowerInfo energies, 1 RawCluster tower-map
// values. sumMembership/momentMembership: 0 the full 7x7 grid, 1 the
// cluster-owned towers. The semantic digest is SHA-256 of the view's semantic
// text (RJ_SHOWER_DEFINITION_FACTORIAL_V1: name, energy source, memberships,
// floor, 7x7 grid, TowerInfo isGood quality, center excluded from the COG-x
// numerator, float32 row-major builder arithmetic); the PhotonID packages bind
// their input view to it.
struct ShowerViewDefinition
{
  const char* name;
  double floorGeV;
  std::int32_t energySource;
  std::int32_t sumMembership;
  std::int32_t momentMembership;
  const char* semanticSha256;
};

inline const std::array<ShowerViewDefinition, 7>& showerViewDefinitions()
{
  static const std::array<ShowerViewDefinition, 7> values{{
      {"H70", 0.070, 0, 0, 1, "fe6633650503f23cf3cfa662a2884d29efbcd7ca465b5b02cfe31ec220db3987"},
      {"H0",  0.0,   0, 0, 1, "bca3665aa23f3f372ab3675d60d2dc9c29973d1da5e0e7053195d563546e2697"},
      {"G70", 0.070, 0, 0, 0, "19be342d81fd569ec9578708e9075fa3192d1baa3a148404855a3411cdc3de77"},
      {"G0",  0.0,   0, 0, 0, "5b45955652de2a98ccc37df1d82593ce3aa7a9e7934861f6bc4ccff1ab5f3a3d"},
      {"O70", 0.070, 0, 1, 1, "044d92f2fe06b495bd97a9215749b33acab92e4ce9553082f3e39fdb670e2755"},
      {"O0",  0.0,   0, 1, 1, "4109b126c999bde784cd09408e6c99160db73f19a8d7a5eaf28026daa7aada1a"},
      {"R70", 0.070, 1, 1, 1, "825959fa8bc8bc97e571cbc48f9f3a534028d80481993e906c320106528dca62"}}};
  return values;
}

inline const ShowerViewDefinition& showerViewDefinition(const std::string& name)
{
  for (const auto& value : showerViewDefinitions())
    if (name == value.name) return value;
  throw std::invalid_argument("unknown shower view definition: " + name);
}

// One named photon model evaluated by the builder during reconstruction, with
// its pinned file and exact input order. The optional domain is the producer's
// invocation range (ET >= minEt, ET <= maxEt, |eta| <= maxAbsEta), not a
// certified physics applicability; a NaN bound leaves that side open.
// unevaluatedScore is written where the model is not evaluated. configuredCut
// is the selection cut the historical recipe applied to this score.
struct PhotonModelConfig
{
  std::string scoreName, file, sha256;
  std::vector<std::string> features;
  float minEt = std::numeric_limits<float>::quiet_NaN();
  float maxEt = std::numeric_limits<float>::quiet_NaN();
  float maxAbsEta = std::numeric_limits<float>::quiet_NaN();
  float unevaluatedScore = std::numeric_limits<float>::quiet_NaN();
  double configuredCut = std::numeric_limits<double>::quiet_NaN();
};

// The producer's photon-identification score, evaluated by the builder on the
// primary shower view and recorded with its historical per-candidate witnesses
// (pho_prodbdt_*). Its model file, digest and input order come from the
// canonical PhotonID package. An empty scoreName means no producer score is
// configured: every candidate then records "model not applicable".
struct ProductionScoreConfig
{
  std::string scoreName;
  std::string modelSha256;
  std::string showerDefinition;

  // Candidates inside [validatedEtMin, validatedEtMax) (and, for Au+Au, a
  // centrality inside [validatedCentralityMin, validatedCentralityMax)) are in
  // the validated domain; finite inputs outside it are diagnostic extrapolation.
  double validatedEtMin = 15.0;
  double validatedEtMax = 35.0;
  double validatedCentralityMin = std::numeric_limits<double>::quiet_NaN();
  double validatedCentralityMax = std::numeric_limits<double>::quiet_NaN();

  // Historical working points of this model, recorded as witnesses inside the
  // validated domain: per photon-ET bin (p+p), or lines in the centrality
  // percent, intercept + slope x centrality (Au+Au). Empty when none apply.
  std::vector<double> workingPointEtEdges;
  std::array<std::vector<double>, 3> workingPoints;  // WP70, WP80, WP90 per ET bin
  bool centralityLinearWorkingPoints = false;
  std::array<double, 3> workingPointIntercepts{};
  std::array<double, 3> workingPointSlopes{};
};

// The active photon preselection recorded per candidate (pho_presel_active),
// in the profile's established variant. reference: e11/e33 < 0.98,
// 0.6 < et1 < 1, 0.8 < e32/e35 < 1 and weta < 0.6. newPPG12: the named NPB
// score above its model's configured cut, and et1, e11/e33 and e32/e35 inside
// their open windows; no weta cut. An empty variant records -1.
struct PreselectionConfig
{
  std::string variant;
  std::string npbScoreName;
  double et1Min = std::numeric_limits<double>::quiet_NaN();
  double et1Max = std::numeric_limits<double>::quiet_NaN();
  double e11OverE33Max = std::numeric_limits<double>::quiet_NaN();
  double e32OverE35Min = std::numeric_limits<double>::quiet_NaN();
  double e32OverE35Max = std::numeric_limits<double>::quiet_NaN();
};

// One isolation working-point line: threshold = interceptGeV + slope * x.
struct IsolationLine
{
  double interceptGeV = std::numeric_limits<double>::quiet_NaN();
  double slope = std::numeric_limits<double>::quiet_NaN();
  double sidebandGapGeV = std::numeric_limits<double>::quiet_NaN();
};

// Isolation witness thresholds recorded with every reconstructed cone.
//
// pp: one photon-ET line per radius. Au+Au: one line in the truncated integer
// centrality per radius, inside [centralityMin, centralityMax); the photon-ET
// fallback line applies otherwise, including a nonfinite centrality. These
// are the historical witnesses; the nominal isolation classes are defined
// separately (PPG12 for pp, the pinned Au+Au freeze for Au+Au).
struct IsolationWitnessConfig
{
  std::vector<IsolationLine> lines;  // aligned with PhotonConfig::isolationRadii
  IsolationLine fallback;
  double centralityMin = std::numeric_limits<double>::quiet_NaN();
  double centralityMax = std::numeric_limits<double>::quiet_NaN();
};

// Nominal photon identification, frozen by the PhotonID package
// (PhotonID/models/nominal.json and each system's nominal manifest). Every
// stored candidate is labelled before writing; labels never change retention.
//
// Au+Au: the model is evaluated on the primary H70 view features inside
// etMin <= ET < etMaxExclusive and centralityMin <= centrality <
// centralityMaxExclusive. pp: the PPG12 producer score (7 GeV minimum, -1
// sentinel) is the nominal score where ET >= etMin; it must agree with a
// recompute of the same model within recomputeUlpLimit float32 ulp without
// changing any region. Regions: tight score > WP70; non-tight 0.1 < score <
// WP80; the band between is excluded; otherwise fail.
struct NominalPhotonIdConfig
{
  std::string packageFile, packageSha256;
  std::string system;  // "auau" or "pp", the package key
  std::string modelFile, modelSha256;
  std::string view;
  std::vector<std::string> featureOrder;
  double etMin = std::numeric_limits<double>::quiet_NaN();
  double etMaxExclusive = std::numeric_limits<double>::quiet_NaN();  // NaN: open
  double centralityMin = std::numeric_limits<double>::quiet_NaN();
  double centralityMaxExclusive = std::numeric_limits<double>::quiet_NaN();
  IsolationLine wp70, wp80;  // Au+Au lines in centrality (sidebandGapGeV unused)
  bool usesProducerScore = false;
  double recomputeUlpLimit = 2.0;
};

// Nominal isolation class definition.
//
// Au+Au: the pinned freeze package; its reconstructed witness at one radius
// and subtraction method against threshold = intercept + slope * centrality,
// with a sideband gap and a centrality domain. pp: PPG12 at R = 0.4 against
// the witness's own stored threshold and sideband (its cone already carries
// the producer's SIM correction). Cone < threshold is isolated, cone above
// the sideband is non-isolated, the band between is excluded, anything
// invalid is -1.
struct NominalIsolationConfig
{
  std::string system;
  std::string freezeFile, freezeSha256, derivationSha256;
  std::string definitionJson;  // pp definition, pinned by its canonical digest
  int witnessIndex = -1;
  double radius = std::numeric_limits<double>::quiet_NaN();
  std::int32_t subtractionMethod = 0;
  IsolationLine threshold;  // Au+Au: centrality line and sideband gap
  double centralityMin = std::numeric_limits<double>::quiet_NaN();
  double centralityMaxExclusive = std::numeric_limits<double>::quiet_NaN();
  bool storedThresholds = false;  // pp
};

// Photon capture.
//
// These bounds are the storage domain, not the measurement selection. The
// analysis window is narrower and lives downstream. Storing the wider domain
// is what lets a working point move without returning to the DST.
struct PhotonConfig
{
  std::vector<PhotonModelConfig> models;
  ProductionScoreConfig productionScore;
  PreselectionConfig preselection;
  IsolationWitnessConfig isolationWitness;
  NominalPhotonIdConfig nominalId;
  NominalIsolationConfig nominalIsolation;
  double minEtGeV = 5.0;
  double maxEtGeV = 40.0;
  double maxAbsEta = 0.7;

  // DATA reconstructed-object vertex domain, cm. The event row is retained
  // outside it; DATA jets additionally require a stored photon. The original
  // SIM retention contract does not apply this DATA vertex cut.
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

  // PPG12 photon-yield source arbitration (pp): the PPG12 tower mask, pinned
  // by digest. An empty file means the arbitration does not apply (Au+Au).
  std::string ppg12TowerMaskFile;
  std::string ppg12TowerMaskSha256;
  std::string ppg12TowerMaskHistogram;
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
  // payloads (the plan's conditions carry them). Embedded simulation consumes
  // the centrality already present in its input instead.
  bool reconstructDataCentrality = false;
  std::string calibrationTag;
};

// Producer-side weight inputs.
//
// Preserve the established producer weight and its factors. PPG12 stores
// relative slice cross section x truth-vertex x SI/DI fraction x period share.
// AuAu vertex and centrality fits are applied downstream without regenerating
// trees; absolute sample normalization also remains downstream.
// Alias factors are never multiplied twice.
struct WeightConfig
{
  bool applyVertexReweight = false;

  std::string vertexReweightFile;
  std::string vertexReweightFileSha256;
  std::string vertexReweightHistogram;
  bool ppg12Nominal = false;
  bool doubleInteraction = false;
  int ppg12LaneCode = 0;  // 1 photon+jet, 2 inclusive jet
  std::string ppg12Slice;  // generator slice of the sample, e.g. photon20 or jet12
  int ppg12SliceCode = 0;  // the slice's bookkeeping sample code
  // Event-ownership window of the slice, [low, high] GeV (high +inf when open).
  double ownershipLowGeV = std::numeric_limits<double>::quiet_NaN();
  double ownershipHighGeV = std::numeric_limits<double>::quiet_NaN();
  // Inclusive-jet slices: photon candidates above this ET belong to the next
  // slice (+inf for photon slices).
  double clusterEtUpperGeV = std::numeric_limits<double>::infinity();
  double sliceFactor = 1.0;
  double interactionFactor = 1.0;
  double periodFactor = 1.0;
};

// One PPG12 crossing-angle period: its run range, luminosity share, SI/DI
// mixture and closure vertex window. An empty key means no period contract;
// the remaining values are then the producer's unset values, which its audit
// rows record as such.
struct PPG12PeriodConfig
{
  std::string key;  // "0mrad" or "1p5mrad"
  int runMin = -1;
  int runMaxExclusive = -1;
  double luminosityPbInv = 1.0;
  double targetLuminosityPbInv = 1.0;
  double doubleFraction = -1.0;
  double closureZCutCm = -1.0;
};

// The PPG12 pp photon-ID lines in photon ET (PhotonID/models/PP/nominal): a
// score above tight() is tight; non-tight lies strictly between
// nonTightLow() and nonTightHigh().
struct PPG12IdLines
{
  static double tight(const double et) { return 0.815625 - 0.0015625 * et; }
  static double nonTightLow(const double et) { return 0.7333333333333333 - 0.01333333333333333 * et; }
  static double nonTightHigh(const double et) { return 0.684375 + 0.0015625 * et; }
};

// The established pp producer's PPG12 bookkeeping (PPG12.cc and
// PPG12Simulation.cc). Photons in these ET bins (GeV, [low, high)) entered its
// per-photon audit. In simulation its histograms also used the recoil-jet
// definition and the unfolding axes below.
struct PPG12BookkeepingConfig
{
  std::vector<double> photonEtBins;

  double recoilJetPtMinGeV = std::numeric_limits<double>::quiet_NaN();
  double recoilDeltaPhiMin = std::numeric_limits<double>::quiet_NaN();  // rad
  double jetMatchMaxDeltaR = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> unfoldRecoPhotonEtBins;
  std::vector<double> unfoldTruthPhotonEtBins;
  std::vector<double> unfoldXJBins;

  // The run the producer seeded its simulation response-ET draws with: its
  // conditions timestamp.
  std::uint64_t responseSeedRun = 0;
};

// The historical pp producer's event classification. It is recorded, never
// applied: it fills evt_terminal_status and the PPG12 vertex/period audit
// occurrences (pdiag_*), and no event or object is removed by it.
//
// DATA: an event is classified when its mandatory MBD vertex exists. The
// configured triggers (GL1 name -> stream key) and the PPG12 scaled-bit audit
// stream define the active streams; then the period run filter and the |vz|
// cut apply. SIM uses the single stream "SIM" and the |vz| cut only.
struct EventClassificationConfig
{
  bool enabled = false;
  std::vector<std::pair<std::string, std::string>> triggers;
  std::string scaledBitStream;
  int scaledBit = -1;
  bool useVertexCut = false;
  double vertexCutCm = std::numeric_limits<double>::quiet_NaN();
  std::vector<PPG12PeriodConfig> periods;  // every period, for run labels
  PPG12PeriodConfig period;                // this job's period (key empty: none)
  bool filterDataByPeriod = false;
  bool useLuminosityWeight = true;
  bool doubleInteraction = false;
  double mixWeight = 1.0;  // SI/DI mixture weight (SIM)
  bool mixWeightFromPeriod = false;
  bool vertexFileFromPeriod = false;

  // Au+Au: the GL1 scaled bits that activate a DATA event (bit, stream), the
  // minimum-bias requirement, the analysis centrality range [low, high) in
  // truncated percent, and each embedded jet sample's ownership window
  // [low, high) GeV of the maximum R = 0.4 truth-jet pT (high +inf when open).
  std::vector<std::pair<int, std::string>> scaledTriggers;
  bool requireMinimumBias = false;
  int centralityMinPercent = 0;
  int centralityMaxPercent = 0;
  std::map<int, std::pair<double, double>> embeddedJetWindows;
};

// The established producer's Au+Au embedded truth-photon misses (truthmiss_*):
// an accepted embedded event whose leading truth signal photon has no
// reconstructed photon with |eta| < absEtaMax, ET in [etMinGeV, etMaxGeV), the
// active preselection, a passing isolation witness at isolationRadius and the
// production's tight identification, or whose leading such photon is another
// truth particle. The tight identification is the established cut-based one:
// weta and wphi in (widthMin, widthMaxIntercept + widthMaxSlopePerGeV x ET),
// e11/e33, et1 and e32/e35 inside their open windows.
struct TruthMissConfig
{
  bool enabled = false;
  double etMinGeV = std::numeric_limits<double>::quiet_NaN();
  double etMaxGeV = std::numeric_limits<double>::quiet_NaN();
  double absEtaMax = std::numeric_limits<double>::quiet_NaN();
  double isolationRadius = std::numeric_limits<double>::quiet_NaN();
  double widthMin = std::numeric_limits<double>::quiet_NaN();
  double widthMaxIntercept = std::numeric_limits<double>::quiet_NaN();
  double widthMaxSlopePerGeV = std::numeric_limits<double>::quiet_NaN();
  std::array<double, 2> e11OverE33{{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()}};
  std::array<double, 2> et1{{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()}};
  std::array<double, 2> e32OverE35{{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()}};
};

// One eventTree branch as the schema declares it.
struct OutputField
{
  std::string name, type, title, counter;
  int length = 0;  // scalar 0, vector -1, fixed array positive
  bool reserved = false;
};

struct OutputConfig
{
  std::vector<OutputField> fields;
  std::string schemaSha256;
  std::vector<std::string> modelFeatures;
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
  std::string resolvedConfigurationSha256;

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
  EventClassificationConfig classification;
  TruthMissConfig truthMisses;
  PPG12BookkeepingConfig ppg12;
  OutputConfig output;
  ProvenanceConfig provenance;

  // Names of the shower-shape definitions to evaluate for every candidate,
  // each one of showerViewDefinitions(). "H70" and "H0" are the hybrid
  // definition (rectangular sums over the full good-tower grid, moments over
  // cluster-owned towers) at a 70 MeV and a 0 MeV floor; these are the two
  // definitions the identification models have used.
  std::vector<std::string> showerDefinitions = {"H70", "H0", "G70", "G0", "O70", "O0", "R70"};

  // Write an event row for every event the producer is handed, including
  // events with no photon candidate. The event population is the exposure
  // denominator of every later normalisation, so it is not filtered by object
  // yield.
  bool retainEveryProducerEncounter = true;
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

// One PPG12 audit occurrence (pdiag_*): a stream and kind with the stage,
// value or vertex the historical bookkeeping recorded, in emission order.
// Unset integers keep their -1/0 defaults and unset measurements stay NaN.
struct PPG12DiagnosticRecord
{
  std::string trigger;
  std::string kind;
  std::int32_t stageCode = 0;
  std::int32_t sourceCode = 0;
  std::int32_t sampleCode = 0;
  std::int32_t decision = 0;
  std::int32_t clusterOrdinal = -1;
  std::int32_t truthTrackId = -1;
  std::int32_t generatorBarcode = -1;
  double sourcePhotonPt = std::numeric_limits<double>::quiet_NaN();
  double flowValue = std::numeric_limits<double>::quiet_NaN();
  double recoPhotonPt = std::numeric_limits<double>::quiet_NaN();
  double responsePhotonPt = std::numeric_limits<double>::quiet_NaN();
  double truthPhotonPt = std::numeric_limits<double>::quiet_NaN();
  double truthPriorWeight = std::numeric_limits<double>::quiet_NaN();
  double recoVertexZ = std::numeric_limits<double>::quiet_NaN();
  double hardTruthVertexZ = std::numeric_limits<double>::quiet_NaN();
  double minimumBiasTruthVertexZ = std::numeric_limits<double>::quiet_NaN();
  double vertexWeight = std::numeric_limits<double>::quiet_NaN();
  double periodEventWeight = std::numeric_limits<double>::quiet_NaN();
  double clusterEnergy = std::numeric_limits<double>::quiet_NaN();
  double clusterEta = std::numeric_limits<double>::quiet_NaN();
  double clusterEt = std::numeric_limits<double>::quiet_NaN();
  double truthEnergy = std::numeric_limits<double>::quiet_NaN();
  double truthEta = std::numeric_limits<double>::quiet_NaN();
  double truthEt = std::numeric_limits<double>::quiet_NaN();
  double energyContribution = std::numeric_limits<double>::quiet_NaN();
  double occurrenceWeight = 1.0;
  std::uint64_t selectionBits = 0;
};

// One Au+Au embedded diagnostic occurrence (embdiag_*): the sample stitch
// decision, an audit category of the truth-isolation diagnostic, or a truth
// photon it accepted, in emission order. Unset fields keep their defaults.
struct EmbeddedDiagnosticRecord
{
  std::string kind;  // STITCH, AUDIT or TRUTH_PHOTON
  std::int32_t auditBin = 0;
  std::int32_t sampleCode = 0;
  std::int32_t centralityBinIndex = -1;
  std::int32_t photonClass = 0;
  double truthPt = std::numeric_limits<double>::quiet_NaN();
  double truthIsolationWitness = std::numeric_limits<double>::quiet_NaN();
  std::int32_t stitchDecision = 0;
  double stitchPhotonPt = std::numeric_limits<double>::quiet_NaN();
};

// One Au+Au embedded truth-photon miss (truthmiss_*): the event's centrality
// class index (0-20, 20-50, 50-80%), the missed truth photon pT, and the
// established producer's decision path (NO_SELECTED_RECO, or
// SELECTED_RECO_NOT_TRUTH_SIGNAL), recorded in the file's constants.
struct TruthMissRecord
{
  std::int32_t centralityBin = -1;
  double truthPt = std::numeric_limits<double>::quiet_NaN();
  std::string decisionPath;
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
  std::uint64_t triggerBunchNumber = std::numeric_limits<std::uint64_t>::max();
  std::uint32_t triggerScalerDiscontinuityFlags = 0;
  std::uint64_t triggerPacketStatus = 0;

  std::int32_t triggerPacketVersion = 0;
  std::uint32_t triggerDecisionsAvailable = 0;
  bool triggerPacketValid = false;
  // Retained GL1 decision contract: 0 valid, 1 missing, 2 unsupported,
  // 3 read error, 4 nonzero status, 5 simulation (not detector data).
  std::int32_t triggerDecisionState = 1;

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

  // Reconstructed vertex satisfies the mode's storage domain. DATA also
  // requires a captured photon before its jet collections are retained.
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
  // Sequential binary64 sum of the finite positive PMT charges; NaN without a
  // PMT container. Complete: 128 channels, each with its identity and charge.
  double mbdTotalCharge = std::numeric_limits<double>::quiet_NaN();
  bool mbdPmtArrayComplete = false;

  bool mbdValid = false;
  bool mbdAvailable = false;
  std::int32_t mbdNativeValid = -1;
  std::int32_t mbdSouthNpmt = -1;
  std::int32_t mbdNorthNpmt = -1;

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

  std::int32_t centralityReplayVersion = 1;
  TriState centralityPmtAvailable = TriState::False;
  std::vector<std::int32_t> centralityPmtId;
  std::vector<double> centralityPmtCharge;
  std::vector<double> centralityPmtTimeNs;
  std::vector<std::int32_t> centralityPmtValid;
  std::vector<std::int32_t> centralityPmtSelected;
  double centralityMbdZ = std::numeric_limits<double>::quiet_NaN();
  bool centralityMbdZValid = false;
  std::int32_t centralityMbdEvent = -1;
  std::int32_t centralityMbdClock = -1;
  std::int32_t centralityMbdFemclock = -1;
  bool centralityInputsValid = false;
  double centralitySelectedCharge = std::numeric_limits<double>::quiet_NaN();

  std::int32_t centralityNativeBin = -1;
  double centralityNativeCentile = std::numeric_limits<double>::quiet_NaN();  // fraction, not percent
  bool centralityNativeValid = false;

  std::int32_t centralityBin = -1;
  // The analysis centrality: the truncated percent of the mbd_NS bin accessor,
  // independent of the coarse jet bin and the native validity witness. pp -1.
  std::int32_t analysisCentrality = -1;
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
  double maximumRawClusterEnergy = std::numeric_limits<double>::quiet_NaN();

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

  // The embedded photon sample's generator-filter witness uses the original
  // HepMC selection, before any reconstructed-photon acceptance. -1 is the
  // measured absence of a qualifying photon; false means capture never ran.
  bool embeddedFilterPhotonCaptured = false;
  double embeddedFilterPhotonPt = -1.0;

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
  // Version/domain stay zero when capture did not enter its vertex domain.
  // Domain 1 is finite photon kinematics, 5 <= ET < 40 and |eta| < 0.7.
  std::int32_t recoPhotonCaptureVersion = 0, recoPhotonCaptureDomain = 0;
  std::uint32_t recoPhotonInspectedCount = 0, recoPhotonEligibleCount = 0;
  std::uint32_t recoPhotonWrittenCount = 0;
  std::uint32_t recoPhotonUnclassifiableCount = 0;

  //
  // PPG12 photon-yield source arbitration (pp with the PPG12 tower mask)
  //
  // State 1 when the arbitration ran over this event's photon clusters, else
  // 0. The leading source is the unmasked cluster with the highest cluster_pt
  // of at least 5 GeV, the first encountered on ties; key 0, encounter -1 and
  // NaN kinematics when there is none.
  //
  std::int32_t ppg12ArbitrationState = 0;
  std::uint32_t ppg12LeadKey = 0;
  std::int32_t ppg12LeadEncounter = -1;
  double ppg12LeadEt = std::numeric_limits<double>::quiet_NaN();
  double ppg12LeadEta = std::numeric_limits<double>::quiet_NaN();
  double ppg12LeadPhi = std::numeric_limits<double>::quiet_NaN();

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

  // Historical event classification (EventClassificationConfig). The
  // mandatory-vertex flag, the active streams in configuration order, the
  // terminal status (pp: 0 none, 1 trigger, 2 vz, 3 period, 4 ownership;
  // Au+Au: 0 none, 1 trigger, 2 vz, 3 minimum bias, 4 stitch, 5 centrality),
  // the PPG12 audit occurrences and the Au+Au embedded diagnostics and truth
  // misses. Diagnostic only; never a retention cut.
  bool nodesReady = false;
  std::vector<std::string> activeTriggers;
  std::int32_t terminalStatus = 0;
  std::vector<PPG12DiagnosticRecord> ppg12Diagnostics;
  std::vector<EmbeddedDiagnosticRecord> embeddedDiagnostics;
  std::vector<TruthMissRecord> truthMisses;

  // pp simulation sample ownership (PPG12 stitching). The stitch inputs are
  // the maximum stored truth photon pT of each established definition, in
  // order g4Stored, hepmcParent, hepmcParentStable, hepmcAnyPhoton, g4NoEmbed,
  // the embedded |eta| < 0.7 source and the maximum truth R=0.4 jet pT; -1
  // when absent. The event weight applies only to an owned, accepted event.
  std::array<double, 5> ppg12StitchPhotonPt{{-1.0, -1.0, -1.0, -1.0, -1.0}};
  double ppg12SourcePhotonPt = -1.0;
  double ppg12MaxTruthJetPtR04 = -1.0;
  bool ppg12Owned = false;
  bool ppg12WeightActive = false;

  // pp simulation: the fills the established producer's weighted and raw
  // histograms received in this event, indexed by weighting lane code (1 to
  // 4; index 0 unused). Fills before the weight configuration are not counted.
  std::array<std::int64_t, 5> ppg12WeightedFills{{0, 0, 0, 0, 0}};
  std::array<std::int64_t, 5> ppg12RawFills{{0, 0, 0, 0, 0}};
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

// Read-only witness of the exact native model invocation. Inputs retain their
// binary32 evaluation representation; object id and container key are distinct.
// Evaluation states follow the retained capture contract: 0 not captured,
// 1 unavailable model, 2 outside domain, 3 invalid input, 4/5 finite/nonfinite
// evaluation. Capturing a witness does not certify model applicability.
struct ModelInvocationRecord
{
  std::int32_t captureVersion = 0, evaluationState = 0, domainMode = 0;
  std::int32_t inProducerDomain = -1, modelAvailable = -1;
  std::int32_t inputOrderIsNpb25 = 0, inputsCaptured = 0, inputsFinite = 0;
  float minEt = std::numeric_limits<float>::quiet_NaN();
  float maxEt = std::numeric_limits<float>::quiet_NaN();
  float maxAbsEta = std::numeric_limits<float>::quiet_NaN();
  float rawScore = std::numeric_limits<float>::quiet_NaN();
  std::uint32_t sourcePhotonId = 0, sourceNativeClusterKey = 0;
  float sourceClusterEt = std::numeric_limits<float>::quiet_NaN();
  float sourceEta = std::numeric_limits<float>::quiet_NaN();
  float sourcePhi = std::numeric_limits<float>::quiet_NaN();
  std::int32_t bindingPolicy = 0, sourceIdentityAvailable = 0;
  std::int32_t kinematicInputBitsAvailable = 0;
  std::uint32_t inputEtBits = 0, inputEtaBits = 0;
  std::int64_t sourceEncounterOrdinal = -1;
  std::string sourceClusterNode, modelFile;
  std::vector<std::string> featureNames;
  std::vector<float> orderedInputs;
  std::vector<std::int32_t> orderedInputFinite;
};

inline const std::vector<std::string>& npbFeatureNames()
{
  static const std::vector<std::string> names={
    "cluster_Et","cluster_Eta","vertexz","e11_over_e33","e32_over_e35",
    "e11_over_e22","e11_over_e13","e11_over_e15","e11_over_e17",
    "e11_over_e31","e11_over_e51","e11_over_e71","e22_over_e33",
    "e22_over_e35","e22_over_e37","e22_over_e53","cluster_weta_cogx",
    "cluster_wphi_cogx","cluster_et1","cluster_et2","cluster_et3",
    "cluster_et4","cluster_w32","cluster_w52","cluster_w72"};
  return names;
}

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

  // Frozen feature-view witnesses; undefined ratios/moments remain NaN.
  double rawCenterEta = std::numeric_limits<double>::quiet_NaN();
  double rawCenterPhi = std::numeric_limits<double>::quiet_NaN();
  double e11OverE33 = std::numeric_limits<double>::quiet_NaN();
  double e32OverE35 = std::numeric_limits<double>::quiet_NaN();
  std::array<double, 6> momentSums{};
  std::int32_t energySource = 0, rectangularMembership = 0, momentMembership = 1;
  std::uint32_t activeSumCellCount = 0, activeMomentCellCount = 0;

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
  // identification model takes these two. The configured photon builder
  // supplies them; an unavailable measurement remains NaN.
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
  bool towerAvailable = false;

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

// The reconstructed isolation witness of one candidate at one radius: the
// cone sum of the system's nominal method (pp topological clusters, Au+Au
// SUB1 towers), its configured working-point threshold and sideband, and
// whether the cone passes. A cone without a valid sum keeps NaN and fails.
struct IsolationWitnessRecord
{
  double radius = std::numeric_limits<double>::quiet_NaN();
  std::int32_t subtractionMethod = 0;   // 1 pp raw, 2 Au+Au SUB1
  std::int32_t reconstructedOrTruth = 0;  // 0 reconstructed
  double coneSum = std::numeric_limits<double>::quiet_NaN();
  double threshold = std::numeric_limits<double>::quiet_NaN();
  double sidebandThreshold = std::numeric_limits<double>::quiet_NaN();
  std::int32_t passState = 0;
  // The historical witness never counted constituents; they are the isoc_ rows.
  std::int32_t constituentCount = 0;
};

// One object entering an isolation cone. The constituents are stored because
// the cone sums cannot otherwise be re-derived with a different convention.
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

  // -1 marks an unmeasured identity, including DATA and a missing primary.
  std::int32_t trackId = -1;
  std::int32_t pid = 0;
  std::int32_t barcode = -1;
  std::int32_t embeddingId = 0;
  std::int32_t vertexId = -1;

  double energyContribution = std::numeric_limits<double>::quiet_NaN();
};

// One cluster of the photon builder's input container as the PPG12
// simulation bookkeeping reads it: its energy; its energy vector from the
// reconstructed vertex (RawClusterUtility::GetEVec), whose magnitude and
// massless pT, eta, phi the producer used; its pseudorapidity from the same
// vertex (GetPseudorapidity); and the Geant4 primary depositing most of its
// energy. Not stored.
struct PPG12ClusterRecord
{
  double energy = std::numeric_limits<double>::quiet_NaN();
  double directionEnergy = std::numeric_limits<double>::quiet_NaN();
  double pt = std::numeric_limits<double>::quiet_NaN();
  double eta = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();
  double vertexEta = std::numeric_limits<double>::quiet_NaN();

  bool primaryFound = false;
  std::int32_t truthTrackId = -1;
  double energyContribution = std::numeric_limits<double>::quiet_NaN();
};

// One reconstructed photon candidate.
//
// Candidates are admitted on kinematics alone. No identification score and no
// isolation working point takes part in retention, which is what lets a model
// change without re-reading the DST.
// The nominal labels of one candidate (NominalPhotonIdConfig and
// NominalIsolationConfig), stored as pho_bdt_score, pho_bdt_state,
// pho_id_class and pho_iso_class.
struct NominalPhotonRecord
{
  double score = std::numeric_limits<double>::quiet_NaN();  // NaN unless scored
  std::int32_t state = 0;           // 1 scored, 0 outside the domain, 2 invalid input
  std::int32_t idClass = -1;        // 3 tight, 2 excluded gap, 1 non-tight, 0 fail, -1 not scored
  std::int32_t isolationClass = -1; // 2 isolated, 1 excluded gap, 0 non-isolated, -1 invalid
};

// Source totals of the nominal labels, written with the file's receipts.
struct NominalReceiptRecord
{
  std::uint64_t photons = 0;
  std::uint64_t parityMismatches = 0;
  std::map<int, std::uint64_t> states, idClasses, isolationClasses;
  double maxRecomputeUlp = 0.0;
  std::uint64_t recomputeRegionFlips = 0;
};

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

  // Copy the builder's primary representation separately from factorial views:
  // its double ratios and a view's float ratios have distinct persisted bits.
  static constexpr std::array<const char*, 15> primaryShapeNames{{
      "weta_cogx", "wphi_cogx", "weta33_cogx", "wphi33_cogx",
      "weta35_cogx", "wphi53_cogx", "et1", "et2", "et3", "et4",
      "e11", "e33", "e32", "e35", "bdt_score"}};
  std::array<double, primaryShapeNames.size()> primaryShapes{};

  // The vertex the reconstruction used to turn cluster position into
  // direction. Stored because it is part of the measurement.
  double producerVertexZ = std::numeric_limits<double>::quiet_NaN();

  // The clusterizer's shower-profile probability (RawCluster::get_prob). It
  // is not a tree branch; the PPG12 bookkeeping's common selection reads it.
  double clusterProbability = std::numeric_limits<double>::quiet_NaN();

  // PPG12 tower mask at this cluster's tower centroid: 1 masked, 0 not
  // masked, -1 when the arbitration does not apply.
  std::int32_t ppg12TowerMaskState = -1;

  PhotonTimingRecord timing;
  DominantTruthRecord dominantTruth;

  ModelInvocationRecord npbCapture, auauNpbCapture;

  // The configured producer score (ProductionScoreConfig) as the builder
  // wrote it on this candidate: the model's unevaluated sentinel outside its
  // invocation domain, NaN when no producer score is configured.
  double productionScore = std::numeric_limits<double>::quiet_NaN();

  std::vector<ShowerShapeRecord> showerShapes;
  std::vector<IsolationRecord> isolation;
  std::vector<IsolationWitnessRecord> isolationWitnesses;

  NominalPhotonRecord nominal;
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
// Original signed constituent energies, accumulated in clustering order.
// Source codes follow Jet::SRC / JetCalo layer semantics, including retowers.
// Widen only the finished float fractions for the final double branches.
struct JetCalorimeterEnergy
{
  std::array<float, 3> energy{{0.F, 0.F, 0.F}};

  static int layer(const int source)
  {
    switch (source)
    {
      case 3: case 4: case 13: case 14: case 17: case 22:
      case 25: case 28: case 29: case 32: case 33: return 0;
      case 5: case 6: case 15: case 18: case 26: case 30:
      case 34: case 35: return 1;
      case 7: case 8: case 16: case 19: case 27: case 31:
      case 36: case 37: return 2;
      default: return -1;
    }
  }

  void add(const int source, const float value)
  {
    const int index = layer(source);
    if (index >= 0) energy[index] += value;
  }

  std::array<double, 3> fractions(const float rawJetEnergy) const
  {
    std::array<double, 3> result;
    result.fill(std::numeric_limits<double>::quiet_NaN());
    if (!std::isfinite(rawJetEnergy) || rawJetEnergy == 0.F) return result;
    for (std::size_t i = 0; i < result.size(); ++i)
      if (std::isfinite(energy[i])) result[i] = static_cast<float>(energy[i] / rawJetEnergy);
    return result;
  }
};

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

  // eta and phi are the pre-JES direction: the raw jet rebuilt through the
  // calibration's binary32 four-vector arithmetic at unit scale, so a JES
  // value never changes the stored direction. calibrationValid is false where
  // the JES function is negative (the calibrated vector is reversed); such a
  // jet keeps its unit-scale pT, and its pairs carry an undefined xJ.
  double rawPt = std::numeric_limits<double>::quiet_NaN();
  double correctedPt = std::numeric_limits<double>::quiet_NaN();
  double eta = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();
  double area = std::numeric_limits<double>::quiet_NaN();

  // EMCal, inner HCal, outer HCal; signed SUB1 values are never clipped.
  std::array<double, 3> calorimeterFractions{{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()}};

  std::uint64_t qualityBitmask = 0;

  bool calibrationValid = false;

  // pT of the unit-scale construction above. Not stored: the established
  // producer's simulation conditions carried no jet energy scale, so its
  // histograms read exactly these (unitScalePt, eta, phi) jets, which the PPG12
  // simulation bookkeeping replays.
  double unitScalePt = std::numeric_limits<double>::quiet_NaN();
};

// One photon and one retained jet of one collection.
//
// The recoil state is the nominal away-side witness, an absolute azimuthal
// separation of at least kRecoilDeltaPhiMin. It is a stored flag and not a
// retention cut: every pair is written, and any away-side requirement is
// reproducible downstream from the stored separation. The rows are the full
// photon x jet product of each event, the largest table in the file.
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
  double energy = std::numeric_limits<double>::quiet_NaN();  // Geant4; not stored

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

  // Raw inventory association is independent of the native signal selection.
  // States follow the final association contract: 0 unknown, 1 no primary,
  // 2 raw photon match, 3 non-photon primary. Truth-only rows carry state 0.
  Identity associationTruthPhotonId;
  std::int32_t associationState = 0;
  double associationDeltaR = std::numeric_limits<double>::quiet_NaN();

  AssociationState state = AssociationState::Unknown;

  double energyContribution = std::numeric_limits<double>::quiet_NaN();
  double deltaR = std::numeric_limits<double>::quiet_NaN();
};

// A reconstructed jet and a truth jet of the same clustering radius.
//
// Every candidate edge within the matching radius is retained, followed by
// separate selected-match rows. Same-radius reconstructed views share one
// truth assignment, preserving the original production matching contract.
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
  std::int64_t physicalEventEnd = -1;

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

  // 0 named bit, 1 TriggerRunInfo node missing, 2 unsupported class,
  // 3 unassigned bit.
  std::int32_t state = 1;

  bool valid = false;
};

}  // namespace photonjet

#endif  // PHOTONJETTREE_INTERNAL_TYPES_H
