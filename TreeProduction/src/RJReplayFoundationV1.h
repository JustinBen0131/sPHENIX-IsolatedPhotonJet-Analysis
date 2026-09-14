#ifndef RJ_REPLAY_FOUNDATION_V1_H
#define RJ_REPLAY_FOUNDATION_V1_H

// Versioned, selection-neutral normalized replay tables for the synchronized
// pp/AuAu photon+jet foundation.  This header is intentionally self-contained
// so the pp and AuAu libraries can still be staged and built independently.

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>
#include <Compression.h>
#include <RtypesCore.h>
#include "RJDominantTruthWitnessV1.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace RJReplayFoundationV1
{
constexpr const char* kSchemaName = "RJ_REPLAY_FOUNDATION_V1";
// Schema14 preserves schema13 raw capture and adds the per-event truth census.
constexpr int kSchemaVersion = 14;
constexpr std::uint64_t kJetFastJetAreaValid = 1ULL << 1U;
constexpr std::uint64_t kJetCalibRawOrdinalLinked = 1ULL << 2U;

// Replay reconstruction must not append into a same-named jet collection
// already carried by an input DST.  A short suffix gives the analysis-owned
// JetReco/JetCalib products a separate node namespace while leaving legacy
// jobs unchanged when the suffix is empty.
inline bool validAnalysisJetNodeSuffix(const std::string& suffix)
{
  if (suffix.empty()) return true;
  if (suffix.size() < 2U || suffix.size() > 32U || suffix.front() != '_')
    return false;
  return std::all_of(suffix.begin() + 1, suffix.end(), [](unsigned char value)
  {
    return std::isalnum(value) != 0 || value == '_';
  });
}

// JetCalib creates one calibrated jet per raw input jet, in input order, and
// sets the calibrated jet id to that zero-based input ordinal.  Raw Jet ids
// are detector/container keys and are not the JetCalib association key.
// Centralize that coresoftware contract so the pp and AuAu writers cannot
// accidentally bind a calibrated jet to a same-numbered but unrelated raw id.
inline bool jetCalibRawOrdinal(std::size_t calibratedSize,
                               std::size_t rawSize,
                               int calibratedOrdinal,
                               int calibratedId,
                               std::size_t& rawOrdinal)
{
  if (calibratedSize != rawSize || calibratedOrdinal < 0 ||
      calibratedId != calibratedOrdinal ||
      static_cast<std::size_t>(calibratedId) >= rawSize)
    return false;
  rawOrdinal = static_cast<std::size_t>(calibratedId);
  return true;
}

// Identity128 remains an internal construction/matching key. New event rows
// use cheap native-scoped identities; makeIdentity remains available for
// immutable dictionary/configuration identities and protected legacy tests.
// Neither representation is persisted as a cryptographic row identity in
// schema 10: the writer emits compact shard-scoped references.
enum class CompactIdentityDomain : std::uint64_t
{
  NULL_REF=0,
  SOURCE=1,
  EVENT=2,
  CANDIDATE=3,
  MODEL=4,
  SHOWER_DEFINITION=5,
  ISOLATION_CONSTITUENT=6,
  ISOLATION_WITNESS=7,
  JET=8,
  JET_CONSTITUENT=9,
  PAIR=10,
  TRUTH_PHOTON=11,
  OCCURRENCE=12,
  TRUTH_JET=13,
  LINK=14,
  SNAPSHOT=15
};

// ROOT leaf-list codes '/l' and '/L' are defined in terms of ULong64_t and
// Long64_t, respectively.  std::uint64_t/std::int64_t are intentionally kept
// for the detector-neutral API and identity arithmetic; on LP64 Linux those
// standard types are not necessarily the same C++ types ROOT expects when a
// branch address is bound.
using SerializedUInt64 = ULong64_t;
using SerializedInt64 = Long64_t;
static_assert(sizeof(SerializedUInt64) == sizeof(std::uint64_t),
              "ROOT unsigned 64-bit persistence type changed width");
static_assert(sizeof(SerializedInt64) == sizeof(std::int64_t),
              "ROOT signed 64-bit persistence type changed width");

struct Identity128
{
  std::uint64_t hi = 0;
  std::uint64_t lo = 0;

  bool isNull() const { return hi == 0 && lo == 0; }
  bool operator==(const Identity128& other) const { return hi == other.hi && lo == other.lo; }
  bool operator!=(const Identity128& other) const { return !(*this == other); }
  std::string hex() const
  {
    std::ostringstream os;
    os << std::hex << std::setfill('0') << std::setw(16) << hi << std::setw(16) << lo;
    return os.str();
  }
};

struct SerializedIdentity128
{
  SerializedUInt64 hi = 0;
  SerializedUInt64 lo = 0;
};

struct IdentityHash
{
  std::size_t operator()(const Identity128& id) const noexcept
  {
    return static_cast<std::size_t>(id.hi ^ (id.lo + 0x9e3779b97f4a7c15ULL + (id.hi << 6U) + (id.hi >> 2U)));
  }
};

namespace detail
{
inline std::uint32_t rotr(std::uint32_t x, std::uint32_t n) { return (x >> n) | (x << (32U - n)); }

inline std::array<std::uint8_t, 32> sha256(const std::string& input)
{
  static constexpr std::uint32_t k[64] = {
    0x428a2f98U,0x71374491U,0xb5c0fbcfU,0xe9b5dba5U,0x3956c25bU,0x59f111f1U,0x923f82a4U,0xab1c5ed5U,
    0xd807aa98U,0x12835b01U,0x243185beU,0x550c7dc3U,0x72be5d74U,0x80deb1feU,0x9bdc06a7U,0xc19bf174U,
    0xe49b69c1U,0xefbe4786U,0x0fc19dc6U,0x240ca1ccU,0x2de92c6fU,0x4a7484aaU,0x5cb0a9dcU,0x76f988daU,
    0x983e5152U,0xa831c66dU,0xb00327c8U,0xbf597fc7U,0xc6e00bf3U,0xd5a79147U,0x06ca6351U,0x14292967U,
    0x27b70a85U,0x2e1b2138U,0x4d2c6dfcU,0x53380d13U,0x650a7354U,0x766a0abbU,0x81c2c92eU,0x92722c85U,
    0xa2bfe8a1U,0xa81a664bU,0xc24b8b70U,0xc76c51a3U,0xd192e819U,0xd6990624U,0xf40e3585U,0x106aa070U,
    0x19a4c116U,0x1e376c08U,0x2748774cU,0x34b0bcb5U,0x391c0cb3U,0x4ed8aa4aU,0x5b9cca4fU,0x682e6ff3U,
    0x748f82eeU,0x78a5636fU,0x84c87814U,0x8cc70208U,0x90befffaU,0xa4506cebU,0xbef9a3f7U,0xc67178f2U};
  std::vector<std::uint8_t> msg(input.begin(), input.end());
  const std::uint64_t bitLen = static_cast<std::uint64_t>(msg.size()) * 8ULL;
  msg.push_back(0x80U);
  while ((msg.size() % 64U) != 56U) msg.push_back(0U);
  for (int i = 7; i >= 0; --i) msg.push_back(static_cast<std::uint8_t>((bitLen >> (i * 8)) & 0xffU));

  std::uint32_t h[8] = {0x6a09e667U,0xbb67ae85U,0x3c6ef372U,0xa54ff53aU,
                        0x510e527fU,0x9b05688cU,0x1f83d9abU,0x5be0cd19U};
  for (std::size_t off = 0; off < msg.size(); off += 64U)
  {
    std::uint32_t w[64] = {};
    for (int i = 0; i < 16; ++i)
      w[i] = (static_cast<std::uint32_t>(msg[off + 4*i]) << 24U) |
             (static_cast<std::uint32_t>(msg[off + 4*i + 1]) << 16U) |
             (static_cast<std::uint32_t>(msg[off + 4*i + 2]) << 8U) |
             static_cast<std::uint32_t>(msg[off + 4*i + 3]);
    for (int i = 16; i < 64; ++i)
    {
      const std::uint32_t s0 = rotr(w[i-15],7U) ^ rotr(w[i-15],18U) ^ (w[i-15] >> 3U);
      const std::uint32_t s1 = rotr(w[i-2],17U) ^ rotr(w[i-2],19U) ^ (w[i-2] >> 10U);
      w[i] = w[i-16] + s0 + w[i-7] + s1;
    }
    std::uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],hh=h[7];
    for (int i = 0; i < 64; ++i)
    {
      const std::uint32_t s1 = rotr(e,6U) ^ rotr(e,11U) ^ rotr(e,25U);
      const std::uint32_t ch = (e & f) ^ ((~e) & g);
      const std::uint32_t t1 = hh + s1 + ch + k[i] + w[i];
      const std::uint32_t s0 = rotr(a,2U) ^ rotr(a,13U) ^ rotr(a,22U);
      const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t t2 = s0 + maj;
      hh=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }
    h[0]+=a; h[1]+=b; h[2]+=c; h[3]+=d; h[4]+=e; h[5]+=f; h[6]+=g; h[7]+=hh;
  }
  std::array<std::uint8_t,32> out{};
  for (int i=0;i<8;++i) for (int j=0;j<4;++j) out[4*i+j]=static_cast<std::uint8_t>((h[i]>>(24-8*j))&0xffU);
  return out;
}
} // namespace detail

inline Identity128 makeIdentity(const std::string& canonicalInput)
{
  const auto digest = detail::sha256(canonicalInput);
  Identity128 id;
  for (int i=0;i<8;++i) id.hi=(id.hi<<8U)|digest[i];
  for (int i=8;i<16;++i) id.lo=(id.lo<<8U)|digest[i];
  if (id.isNull()) id.lo=1; // reserved null identity; SHA-256 collision guard.
  return id;
}

inline Identity128 makeScopedIdentity(
    CompactIdentityDomain domain,std::uint64_t scope,std::uint64_t ordinal)
{
  constexpr std::uint64_t kScopeMask=(1ULL<<56U)-1ULL;
  const auto code=static_cast<std::uint64_t>(domain);
  if(code==0||code>0xffULL||scope>kScopeMask||
     ordinal==std::numeric_limits<std::uint64_t>::max())
    throw std::invalid_argument("compact transient identity is out of range");
  return Identity128{(code<<56U)|scope,ordinal+1ULL};
}

inline Identity128 identityFromSha256(const std::string& value)
{
  if(value.size()!=64)throw std::invalid_argument("identity SHA-256 must contain 64 hexadecimal characters");
  auto nibble=[](char c)->std::uint64_t
  {
    if(c>='0'&&c<='9')return static_cast<std::uint64_t>(c-'0');
    if(c>='a'&&c<='f')return static_cast<std::uint64_t>(c-'a'+10);
    throw std::invalid_argument("identity SHA-256 must be lowercase hexadecimal");
  };
  Identity128 id;
  for(std::size_t i=0;i<16;++i)
  {
    const auto valueByte=(nibble(value[2*i])<<4U)|nibble(value[2*i+1]);
    if(i<8)id.hi=(id.hi<<8U)|valueByte;
    else id.lo=(id.lo<<8U)|valueByte;
  }
  if(id.isNull())id.lo=1;
  return id;
}

inline std::string sha256Hex(const std::string& canonicalInput)
{
  const auto digest = detail::sha256(canonicalInput);
  std::ostringstream os;
  os << std::hex << std::setfill('0');
  for (const auto byte : digest) os << std::setw(2) << static_cast<unsigned int>(byte);
  return os.str();
}

enum class ModelApplicability : std::int32_t { VALIDATED_DOMAIN=0, DIAGNOSTIC_EXTRAPOLATION=1, MODEL_NOT_APPLICABLE=2, INPUT_INVALID=3 };
enum class LinkClass : std::int32_t { MATCH=0, RECO_FAKE=1, TRUTH_MISS=2, WRONG_PHOTON=3, WRONG_RECOIL=4, MATCH_CANDIDATE=5 };
enum class RecoTruthType : std::int32_t { NONE=0, PHOTON=1, JET=2 };
enum class WeightTargetType : std::int32_t { EVENT=1, CANDIDATE=2, JET=3, TRUTH_PHOTON=4, TRUTH_JET=5 };

struct SourceOccurrenceRow { Identity128 id; std::string lane,dataset,sample,period,si_di_role,ownership_state,input_uri_hash,input_file_sha256,source_manifest_sha256; std::int32_t run=0,segment=0; };

// Selection-neutral accounting for the one native embedded-primary photon
// traversal. A valid empty source is capture_state==1 with all counts zero;
// missing truth input is capture_state==0 and can never certify a zero signal.
struct TruthPhotonInventoryV1
{
  // -2 not recorded, -1 DATA/not applicable, 0 missing/invalid truth input,
  // 1 captured from a valid truth container.
  std::int32_t capture_state=-2;
  std::int32_t native_embedded_primary_pid22_count=0;
  std::int32_t serialized_raw_count=0;
  std::int32_t rejected_count=0;
  std::int32_t rejected_invalid_kinematics_count=0;
  std::int32_t rejected_nonfinite_isolation_count=0;
  std::int32_t duplicate_track_count=0;
  // Counts unusable constituent contributions encountered while evaluating
  // photon isolation. These do not change the nominal photon population, but
  // any nonzero value keeps the denominator explicitly incomplete.
  std::int32_t isolation_input_incomplete_count=0;
  std::int32_t analysis_signal_count=0;

  void rejectInvalidKinematics()
  { ++rejected_count; ++rejected_invalid_kinematics_count; }
  void rejectNonfiniteIsolation()
  { ++rejected_count; ++rejected_nonfinite_isolation_count; }
  void rejectDuplicateTrack()
  { ++rejected_count; ++duplicate_track_count; }
  bool rejectedPartitionCloses() const
  {
    return rejected_count>=0 && rejected_invalid_kinematics_count>=0 &&
        rejected_nonfinite_isolation_count>=0 && duplicate_track_count>=0 &&
        rejected_count==rejected_invalid_kinematics_count+
        rejected_nonfinite_isolation_count+duplicate_track_count;
  }
  bool rawEquationCloses() const
  {
    return native_embedded_primary_pid22_count>=0 && serialized_raw_count>=0 &&
        analysis_signal_count>=0 && analysis_signal_count<=serialized_raw_count &&
        rejectedPartitionCloses() &&
        native_embedded_primary_pid22_count==serialized_raw_count+rejected_count;
  }
  bool denominatorComplete() const
  {
    return capture_state==1 && rawEquationCloses() &&
        isolation_input_incomplete_count==0 && duplicate_track_count==0 &&
        rejected_nonfinite_isolation_count==0;
  }
};

struct EventRow
{
  Identity128 id,source_id;
  std::int32_t run=0;
  std::int64_t event_sequence=0;
  std::int64_t source_file_ordinal=-1,source_entry_ordinal=-1,source_global_entry_ordinal=-1;
  std::int64_t physical_event_sequence=-1;
  std::int32_t physical_event_sequence_valid=0;
  std::uint64_t trigger_bits=0,live_trigger_bits=0,scaled_trigger_bits=0;
  std::uint64_t trigger_packet_status=0;
  int trigger_capture_version=0,trigger_decision_state=1,trigger_packet_version=0;
  unsigned trigger_decision_available=0;
  // Source-local reference to a lossless cumulative-counter snapshot. Zero
  // means not applicable (SIM); DATA missing packets have explicit state rows.
  std::uint64_t trigger_scaler_snapshot_id=0;
  double vertex_z=0,centrality=-1,event_weight=1;
  // Raw event observables needed by already-existing direct scientific
  // objects.  They extend RJEventV1 without adding a seventeenth production
  // tree or changing any selection, weight, histogram, or reconstruction.
  double truth_vertex_z=std::numeric_limits<double>::quiet_NaN();
  double truth_mb_vertex_z=std::numeric_limits<double>::quiet_NaN();
  // -2 not recorded, -1 DATA/not applicable, 0 missing container, 1 captured.
  // Native IDs/embedding flags are evidence, not inferred hard/MB role labels.
  std::int32_t truth_vertex_capture_state=-2,truth_primary_vertex_id=-1;
  std::vector<int> truth_vertex_id,truth_vertex_embedding_id;
  std::vector<int> truth_vertex_embedding_valid,truth_vertex_z_valid;
  std::vector<double> truth_vertex_z_values;
  std::int32_t truth_jet_capture_state=-2;
  std::vector<int> truth_jet_radius_code,truth_jet_container_valid;
  // Event-specific executed views, including empty reconstructed containers.
  // An absent entry is unknown/not captured, never a proven zero-jet view.
  std::vector<int> reco_jet_view_radius_code,reco_jet_view_available;
  std::vector<std::string> reco_jet_view_input_identity,reco_jet_view_subtraction_identity;
  double mbd_total_charge=std::numeric_limits<double>::quiet_NaN();
  std::int32_t centrality_replay_version=-1;
  std::int32_t centrality_pmt_available=-1;
  std::int32_t centrality_inputs_valid=-1;
  std::int32_t centrality_mb_decision=-1;
  std::int32_t centrality_mbd_z_valid=-1;
  std::int32_t centrality_mbd_event=-1;
  std::int32_t centrality_mbd_clock=-1;
  std::int32_t centrality_mbd_femclock=-1;
  std::int32_t centrality_native_bin=-1;
  std::int32_t centrality_native_valid=-1;
  double centrality_mbd_z=std::numeric_limits<double>::quiet_NaN();
  double centrality_selected_charge=std::numeric_limits<double>::quiet_NaN();
  double centrality_native_centile=std::numeric_limits<double>::quiet_NaN();
  std::vector<int> centrality_pmt_id;
  std::vector<double> centrality_pmt_charge;
  std::vector<double> centrality_pmt_time;
  std::vector<int> centrality_pmt_valid;
  std::vector<int> centrality_pmt_selected;
  std::int32_t collaborator_interface_version=2;
  std::int32_t reco_vertex_valid=-1,sample_weight_valid=-1;
  std::int32_t vertex_weight_valid=-1,truth_denominator_complete=-1;
  TruthPhotonInventoryV1 truth_photon_inventory;
  std::int32_t mbd_pmt_available=-1;
  double legacy_event_weight=std::numeric_limits<double>::quiet_NaN();
  std::vector<int> mbd_pmt_id,mbd_pmt_arm,mbd_pmt_charge_valid;
  std::vector<double> mbd_pmt_charge;
  double emcal_total_energy=std::numeric_limits<double>::quiet_NaN();
  double ihcal_total_energy=std::numeric_limits<double>::quiet_NaN();
  double ohcal_total_energy=std::numeric_limits<double>::quiet_NaN();
  double total_calo_energy=std::numeric_limits<double>::quiet_NaN();
  double max_cluster_energy=std::numeric_limits<double>::quiet_NaN();
  std::int32_t event_calo_require_isgood=-1,embedded_minbias_decision=-1;
  // Frozen PPG12 leading-source-cluster arbitration.  This scan occurs after
  // the nominal tower mask and ET >= 5 GeV gate but before the candidate eta
  // cut, exactly as in RecoEffCalculator_TTreeReader.C.  A higher-ET cluster
  // outside |eta|<0.7 can therefore be represented without serializing it as
  // an analysis candidate.
  std::int32_t ppg12_source_arbitration_state=0;
  std::uint32_t ppg12_leading_source_cluster_key=0;
  std::int32_t ppg12_leading_source_encounter_ordinal=-1;
  double ppg12_leading_source_et=std::numeric_limits<double>::quiet_NaN();
  double ppg12_leading_source_eta=std::numeric_limits<double>::quiet_NaN();
  double ppg12_leading_source_phi=std::numeric_limits<double>::quiet_NaN();
  // Compact, selection-neutral witness for the Au+Au background state used
  // when the reconstructed SUB1 jets were made.  p+p rows retain the explicit
  // not-applicable defaults.  The integer contract state avoids repeating a
  // constant node-name string in every event: 0 means not applicable (p+p),
  // 1 means TowerInfoBackground_Sub2 (Au+Au), as sealed by campaign metadata.
  // The UE vectors are the three calorimeter-layer
  // eta profiles stored by TowerInfoBackground_Sub2; no tower payload is
  // duplicated into the replay tables.
  std::int32_t jet_background_contract_state=0;
  std::int32_t jet_background_flow_mode=-1,jet_background_eta_bins=0;
  double jet_background_v2=std::numeric_limits<double>::quiet_NaN();
  double jet_background_psi2=std::numeric_limits<double>::quiet_NaN();
  std::vector<float> jet_background_ue_emcal;
  std::vector<float> jet_background_ue_ihcal;
  std::vector<float> jet_background_ue_ohcal;
  std::int32_t jet_background_n_strips=-1,jet_background_n_towers=-1;
  std::int32_t jet_background_flow_failure_state=-1;
  std::int32_t nodes_ready=0;
  std::int32_t ppg12_weight_active_state=0,ppg12_weight_lane_component_code=0;
  std::int64_t ppg12_weighted_fill_count=0,ppg12_raw_fill_count=0;
  std::array<std::int64_t,5> ppg12_weighted_fill_count_by_code{};
  std::array<std::int64_t,5> ppg12_raw_fill_count_by_code{};
  std::int32_t terminal_status=0,candidate_count=0,tag_count=0,recoil_count=0;
};
// A worker-local encounter counter is not a physical event number. Capture
// the current node explicitly; missing headers remain unknown, never zero or
// a silently substituted shard-local ordinal.
template<class Header>
inline void capturePhysicalEventSequence(EventRow& event, const Header* header)
{
  event.physical_event_sequence=header
      ?static_cast<std::int64_t>(header->get_EvtSequence()):-1;
  event.physical_event_sequence_valid=
      header && event.physical_event_sequence>=0 ? 1 : 0;
}

template<class TruthInfo>
inline void captureTruthVertices(EventRow& event, bool simulation, TruthInfo* truth)
{
  event.truth_vertex_capture_state=simulation ? (truth ? 1 : 0) : -1;
  event.truth_jet_capture_state=simulation ? 0 : -1;
  event.truth_jet_radius_code.clear();event.truth_jet_container_valid.clear();
  event.truth_primary_vertex_id=-1;
  event.truth_vertex_id.clear();event.truth_vertex_embedding_id.clear();
  event.truth_vertex_embedding_valid.clear();event.truth_vertex_z_valid.clear();
  event.truth_vertex_z_values.clear();
  if(!simulation || !truth)return;
  event.truth_primary_vertex_id=truth->GetPrimaryVertexIndex();
  const auto vertices=truth->GetPrimaryVtxRange();
  const auto embeddings=truth->GetEmbeddedVtxIds();
  for(auto it=vertices.first;it!=vertices.second;++it)
  {
    const int id=it->first;
    int embedding=0,embeddingValid=0;
    for(auto e=embeddings.first;e!=embeddings.second;++e)
      if(e->first==id){embedding=e->second;embeddingValid=1;break;}
    const double z=it->second ? it->second->get_z()
        :std::numeric_limits<double>::quiet_NaN();
    event.truth_vertex_id.push_back(id);
    event.truth_vertex_embedding_id.push_back(embedding);
    event.truth_vertex_embedding_valid.push_back(embeddingValid);
    event.truth_vertex_z_values.push_back(z);
    event.truth_vertex_z_valid.push_back(it->second && std::isfinite(z) ? 1 : 0);
  }
}

// Current-event truth nodes, independent of reconstructed-object acceptance.
// Preserve fetchNodes candidate order without cached previous-event pointers.
template<class Lookup>
inline auto captureTruthJetContainers(const std::vector<std::string>& keys,
                                     std::string overrideNode, Lookup lookup)
{
  while(!overrideNode.empty() && std::isspace(static_cast<unsigned char>(overrideNode.front())))
    overrideNode.erase(overrideNode.begin());
  while(!overrideNode.empty() && std::isspace(static_cast<unsigned char>(overrideNode.back())))
    overrideNode.pop_back();
  std::map<std::string,decltype(lookup(std::string{}))> result;
  for(const auto& key:keys)
  {
    if(key.empty() || result.count(key))continue;
    std::vector<std::string> nodes;
    if(!overrideNode.empty())
    {
      const auto placeholder=overrideNode.find("{rKey}");
      if(placeholder!=std::string::npos)
      {
        auto expanded=overrideNode;expanded.replace(placeholder,6,key);
        nodes.push_back(expanded);
      }
      else
      {
        nodes.push_back(overrideNode);
        if(overrideNode.find("_"+key)==std::string::npos)
          nodes.push_back(overrideNode+"_"+key);
      }
    }
    nodes.push_back("AntiKt_Truth_"+key);
    nodes.push_back("AntiKt_TruthFromParticles_"+key);
    result[key]=nullptr;
    for(const auto& node:nodes)
      if(auto* found=lookup(node)){result[key]=found;break;}
  }
  return result;
}

template<class Containers>
inline void recordTruthJetAvailability(EventRow& event,const Containers& containers)
{
  event.truth_jet_capture_state=1;
  event.truth_jet_radius_code.clear();event.truth_jet_container_valid.clear();
  for(const auto& item:containers)
  {
    event.truth_jet_radius_code.push_back(std::stoi(item.first.substr(1)));
    event.truth_jet_container_valid.push_back(item.second ? 1 : 0);
  }
}

struct PhotonCandidateRow
{
  Identity128 id,event_id;
  // Native PhotonClusterContainer/RawCluster key within this event.
  std::uint32_t native_cluster_key=0;
  std::int32_t encounter_ordinal=0,finite_feature_state=0,below15_retention_state=0;
  std::int32_t reference_preselection_state=0,active_preselection_state=0;
  std::int32_t ppg12_tower_mask_state=-1,ppg12_source_eligible_state=0;
  std::int32_t truth_signal_match_state=-1,truth_signal_match_barcode=-1;
  RJDominantTruthWitnessV1::Witness dominant_truth;
  std::vector<double> rank_keys;
  double cluster_et=0,eta=0,phi=0;
  double native_weta_cogx=std::numeric_limits<double>::quiet_NaN();
  double native_wphi_cogx=std::numeric_limits<double>::quiet_NaN();
  double native_weta33_cogx=std::numeric_limits<double>::quiet_NaN();
  double native_wphi33_cogx=std::numeric_limits<double>::quiet_NaN();
  double native_weta35_cogx=std::numeric_limits<double>::quiet_NaN();
  double native_wphi53_cogx=std::numeric_limits<double>::quiet_NaN();
  double native_et1=std::numeric_limits<double>::quiet_NaN();
  double native_e11_over_e33=std::numeric_limits<double>::quiet_NaN();
  double native_e32_over_e35=std::numeric_limits<double>::quiet_NaN();
  double npb_score=std::numeric_limits<double>::quiet_NaN();
  double auau_npb_score=std::numeric_limits<double>::quiet_NaN();
  std::vector<float> ordered_features;
  std::vector<std::string> shower_definition_views;
  std::uint64_t preselection_bitmask=0;
};
struct ModelEvaluationRow { Identity128 candidate_id,model_id; std::string model_sha256,shower_definition_id,shower_semantic_sha256; std::vector<float> ordered_input_witnesses; double raw_score=std::numeric_limits<double>::quiet_NaN(); std::int32_t finite_score=0,applicability_state=2; double wp70=std::numeric_limits<double>::quiet_NaN(),wp80=std::numeric_limits<double>::quiet_NaN(),wp90=std::numeric_limits<double>::quiet_NaN(),delta_wp70=std::numeric_limits<double>::quiet_NaN(),delta_wp80=std::numeric_limits<double>::quiet_NaN(),delta_wp90=std::numeric_limits<double>::quiet_NaN(); };
struct ShowerCellRow { Identity128 candidate_id; std::int32_t local_eta_index=0,local_phi_index=0,tower_eta_index=-1,tower_phi_index=-1; std::uint64_t tower_key=0; double calibrated_energy=0,rawcluster_map_value=std::numeric_limits<double>::quiet_NaN(); std::int32_t is_good=0,is_zero=0,is_negative=0,is_nonfinite=0,seed_state=0,denominator_membership=0,rawcluster_owned=0,rawcluster_value_present=0,floor0_membership=0,floor70_membership=0,grid_membership_bitmask=0; };
struct ShowerFeatureViewRow
{
  Identity128 candidate_id,definition_id;
  std::string definition_name,semantic_sha256;
  std::vector<float> ordered_features;
  double floor_gev=0,cog_eta=std::numeric_limits<double>::quiet_NaN(),cog_phi=std::numeric_limits<double>::quiet_NaN();
  double raw_center_eta=std::numeric_limits<double>::quiet_NaN(),raw_center_phi=std::numeric_limits<double>::quiet_NaN();
  std::int32_t center_eta_index=-1,center_phi_index=-1;
  double e11=0,e33=0,e32=0,e35=0,e11_over_e33=std::numeric_limits<double>::quiet_NaN(),e32_over_e35=std::numeric_limits<double>::quiet_NaN();
  double weta_cogx=std::numeric_limits<double>::quiet_NaN(),wphi_cogx=std::numeric_limits<double>::quiet_NaN(),weta33_cogx=std::numeric_limits<double>::quiet_NaN(),wphi33_cogx=std::numeric_limits<double>::quiet_NaN();
  double native_et1=std::numeric_limits<double>::quiet_NaN(),native_et2=std::numeric_limits<double>::quiet_NaN(),native_et3=std::numeric_limits<double>::quiet_NaN(),native_et4=std::numeric_limits<double>::quiet_NaN();
  double moment_eta_numerator=0,moment_phi_numerator=0,moment_denominator=0,moment33_eta_numerator=0,moment33_phi_numerator=0,moment33_denominator=0;
  std::int32_t energy_source=0,rectangular_membership=0,moment_membership=0,finite_feature_state=0;
  std::int32_t good_cell_count=0,owned_cell_count=0,active_sum_cell_count=0,active_moment_cell_count=0,exact_zero_count=0,negative_count=0,nonfinite_count=0;
};
struct IsolationConstituentRow { Identity128 candidate_id; std::int32_t native_source=-1; std::uint32_t native_key=0; double delta_eta=0,delta_phi=0,delta_r=0,raw_energy=0,calibrated_energy=0,sub1_energy=0,phosub_residual=0; std::int32_t subsystem=0,quality_state=0,mask_state=0,candidate_removal_state=0; };
struct IsolationWitnessRow { Identity128 candidate_id,isolation_id; double radius=0,cone_sum=0,threshold=0,sideband_threshold=0; std::int32_t subtraction_method=0,reconstructed_or_truth=0,pass_state=0,constituent_count=0; };
struct JetRow { Identity128 id,event_id; std::uint32_t native_jet_key=0,native_raw_jet_key=0; std::string algorithm,input_identity,subtraction_identity; double radius=0,raw_pt=0,corrected_pt=0,eta=0,phi=0,mass=std::numeric_limits<double>::quiet_NaN(),area=std::numeric_limits<double>::quiet_NaN(); std::uint64_t quality_bitmask=0; std::int32_t deterministic_order=0; };
struct JetConstituentRow { Identity128 jet_id; std::int32_t native_source=-1; std::uint32_t native_key=0; std::int32_t constituent_ordinal=0,subsystem=0,quality_state=0; double energy=0,eta=0,phi=0; };
struct PhotonJetPairRow { Identity128 id,event_id,candidate_id,jet_id; double delta_phi=0,xjgamma=0; std::int32_t recoil_state=0,photon_rank=0,jet_rank=0,wrong_photon_class=0,wrong_recoil_class=0; };
// source_role retains its legacy analysis-signal indicator (1 or 0).
// sample_source_role independently records photon/jet family role (1/2/unknown 0).
struct TruthPhotonRow { std::int32_t generator_occurrence_embedding_id=-1,sample_source_role=0; Identity128 id,event_id; double pt=0,eta=0,phi=0,truth_isolation_witness=0,truth_isolation_r03=0,truth_isolation_r04=0; std::int32_t prompt_class=0,source_role=0,g4_photon_valid=0,hepmc_association_valid=0,truth_isolation_valid=0,analysis_signal_r03=0,generator_barcode=-1,native_track_id=-1,native_vertex_id=-1,embedding_id=0,reporting_guard_state=0; };
struct TruthPhotonMissOccurrenceRow
{
  Identity128 id,event_id;
  std::string trigger,iso_view,decision_path;
  double truth_pt=0;
  std::int32_t centrality_bin_index=-1,occurrence_ordinal=0;
};
struct EmbeddedPhotonDiagnosticOccurrenceRow
{
  Identity128 id,event_id;
  std::string trigger,occurrence_kind;
  double truth_pt=std::numeric_limits<double>::quiet_NaN();
  double truth_isolation_witness=std::numeric_limits<double>::quiet_NaN();
  double stitch_photon_pt=std::numeric_limits<double>::quiet_NaN();
  std::int32_t audit_bin=0,sample_code=0,centrality_bin_index=-1;
  std::int32_t photon_class=0,stitch_decision=0,occurrence_ordinal=0;
};
struct PPG12DiagnosticOccurrenceRow
{
  Identity128 id,event_id;
  std::string trigger,occurrence_kind;
  std::int32_t stage_code=0,source_code=0,sample_code=0,decision=0;
  std::int32_t cluster_ordinal=-1,truth_track_id=-1,generator_barcode=-1;
  double source_photon_pt=std::numeric_limits<double>::quiet_NaN();
  double flow_value=std::numeric_limits<double>::quiet_NaN();
  double reco_photon_pt=std::numeric_limits<double>::quiet_NaN();
  double response_photon_pt=std::numeric_limits<double>::quiet_NaN();
  double truth_photon_pt=std::numeric_limits<double>::quiet_NaN();
  double truth_prior_weight=std::numeric_limits<double>::quiet_NaN();
  double reco_vertex_z=std::numeric_limits<double>::quiet_NaN();
  double hard_truth_vertex_z=std::numeric_limits<double>::quiet_NaN();
  double mb_truth_vertex_z=std::numeric_limits<double>::quiet_NaN();
  double vertex_weight=std::numeric_limits<double>::quiet_NaN();
  double period_event_weight=std::numeric_limits<double>::quiet_NaN();
  double cluster_energy=std::numeric_limits<double>::quiet_NaN();
  double cluster_eta=std::numeric_limits<double>::quiet_NaN();
  double cluster_et=std::numeric_limits<double>::quiet_NaN();
  double truth_energy=std::numeric_limits<double>::quiet_NaN();
  double truth_eta=std::numeric_limits<double>::quiet_NaN();
  double truth_et=std::numeric_limits<double>::quiet_NaN();
  double energy_contribution=std::numeric_limits<double>::quiet_NaN();
  double occurrence_weight=1.0;
  std::uint64_t selection_bitmask=0;
  std::int32_t occurrence_ordinal=0;
};
struct TruthJetRow { Identity128 id,event_id; std::uint32_t native_jet_key=0; std::string algorithm,ownership_state; double radius=0,pt=0,eta=0,phi=0; std::int32_t reporting_guard_state=0; };
struct RecoTruthLinkRow { Identity128 id,reco_id,truth_id; std::int32_t reco_type=0,truth_type=0,link_class=0; double match_metric=std::numeric_limits<double>::quiet_NaN(); };
struct WeightComponentRow { Identity128 target_id; std::string component_type; double slice_weight=1,cross_section_weight=1,vertex_weight=1,si_di_weight=1,period_weight=1,exposure_weight=1,final_weight=1; std::int32_t target_type=static_cast<std::int32_t>(WeightTargetType::EVENT),application_count=0; };
struct EventDisplaySnapshotRow { Identity128 id,event_id; std::string selection_reason,quota_class,serialized_payload_hash; };

// Selection-neutral jet matching shared by the p+p and Au+Au replay writers.
// Every same-radius edge inside the frozen candidate gate is retained as a
// MATCH_CANDIDATE row.  MATCH/RECO_FAKE/TRUTH_MISS rows then form an exact,
// disjoint final partition.  Candidate rows are diagnostics and are excluded
// from that partition.
struct JetMatchObject
{
  Identity128 id;
  double pt=0,eta=0,phi=0,radius=0;
};

inline std::vector<RecoTruthLinkRow> buildDeterministicJetLinks(
    const std::vector<JetMatchObject>& reco,
    const std::vector<JetMatchObject>& truth,
    double deltaRMax,
    std::uint64_t identityScope=0,
    std::uint64_t firstLinkOrdinal=0)
{
  if (!std::isfinite(deltaRMax) || deltaRMax <= 0.0)
    throw std::invalid_argument("jet match deltaR must be finite and positive");

  struct Candidate
  {
    std::size_t recoIndex=0,truthIndex=0;
    double deltaR=0;
  };
  std::vector<Candidate> candidates;
  std::vector<double> nearestReco(reco.size(),std::numeric_limits<double>::quiet_NaN());
  std::vector<double> nearestTruth(truth.size(),std::numeric_limits<double>::quiet_NaN());
  auto identityLess=[](const Identity128& a,const Identity128& b)
  { return a.hi<b.hi || (a.hi==b.hi && a.lo<b.lo); };
  auto updateNearest=[](double& current,double value)
  { if(!std::isfinite(current)||value<current) current=value; };
  auto wrappedDeltaPhi=[](double a,double b)
  { return std::atan2(std::sin(a-b),std::cos(a-b)); };

  for(std::size_t ir=0;ir<reco.size();++ir)
  {
    if(reco[ir].id.isNull()||!std::isfinite(reco[ir].eta)||!std::isfinite(reco[ir].phi))continue;
    for(std::size_t it=0;it<truth.size();++it)
    {
      if(truth[it].id.isNull()||!std::isfinite(truth[it].eta)||!std::isfinite(truth[it].phi))continue;
      if(std::fabs(reco[ir].radius-truth[it].radius)>1.0e-6)continue;
      const double dr=std::hypot(reco[ir].eta-truth[it].eta,
                                 wrappedDeltaPhi(reco[ir].phi,truth[it].phi));
      if(!std::isfinite(dr))continue;
      updateNearest(nearestReco[ir],dr);
      updateNearest(nearestTruth[it],dr);
      if(dr<deltaRMax)candidates.push_back({ir,it,dr});
    }
  }
  std::sort(candidates.begin(),candidates.end(),[&](const Candidate& a,const Candidate& b)
  {
    if(a.deltaR!=b.deltaR)return a.deltaR<b.deltaR;
    if(reco[a.recoIndex].pt!=reco[b.recoIndex].pt)
      return reco[a.recoIndex].pt>reco[b.recoIndex].pt;
    if(reco[a.recoIndex].id!=reco[b.recoIndex].id)
      return identityLess(reco[a.recoIndex].id,reco[b.recoIndex].id);
    return identityLess(truth[a.truthIndex].id,truth[b.truthIndex].id);
  });

  std::vector<RecoTruthLinkRow> links;
  links.reserve(candidates.size()+reco.size()+truth.size());
  for(const Candidate& edge:candidates)
  {
    RecoTruthLinkRow row;
    row.id=makeScopedIdentity(CompactIdentityDomain::LINK,identityScope,
                              firstLinkOrdinal+links.size());
    row.reco_type=static_cast<int>(RecoTruthType::JET);
    row.reco_id=reco[edge.recoIndex].id;
    row.truth_type=static_cast<int>(RecoTruthType::JET);
    row.truth_id=truth[edge.truthIndex].id;
    row.match_metric=edge.deltaR;
    row.link_class=static_cast<int>(LinkClass::MATCH_CANDIDATE);
    links.push_back(row);
  }

  std::unordered_set<Identity128,IdentityHash> matchedReco,matchedTruth;
  for(const Candidate& edge:candidates)
  {
    const auto& recoObject=reco[edge.recoIndex];
    const auto& truthObject=truth[edge.truthIndex];
    if(matchedReco.count(recoObject.id)||matchedTruth.count(truthObject.id))continue;
    matchedReco.insert(recoObject.id);matchedTruth.insert(truthObject.id);
    RecoTruthLinkRow row;
    row.id=makeScopedIdentity(CompactIdentityDomain::LINK,identityScope,
                              firstLinkOrdinal+links.size());
    row.reco_type=static_cast<int>(RecoTruthType::JET);row.reco_id=recoObject.id;
    row.truth_type=static_cast<int>(RecoTruthType::JET);row.truth_id=truthObject.id;
    row.match_metric=edge.deltaR;row.link_class=static_cast<int>(LinkClass::MATCH);
    links.push_back(row);
  }
  for(std::size_t ir=0;ir<reco.size();++ir)
  {
    if(matchedReco.count(reco[ir].id))continue;
    RecoTruthLinkRow row;
    row.id=makeScopedIdentity(CompactIdentityDomain::LINK,identityScope,
                              firstLinkOrdinal+links.size());
    row.reco_type=static_cast<int>(RecoTruthType::JET);row.reco_id=reco[ir].id;
    row.truth_type=static_cast<int>(RecoTruthType::NONE);
    row.match_metric=nearestReco[ir];row.link_class=static_cast<int>(LinkClass::RECO_FAKE);
    links.push_back(row);
  }
  for(std::size_t it=0;it<truth.size();++it)
  {
    if(matchedTruth.count(truth[it].id))continue;
    RecoTruthLinkRow row;
    row.id=makeScopedIdentity(CompactIdentityDomain::LINK,identityScope,
                              firstLinkOrdinal+links.size());
    row.reco_type=static_cast<int>(RecoTruthType::NONE);
    row.truth_type=static_cast<int>(RecoTruthType::JET);row.truth_id=truth[it].id;
    row.match_metric=nearestTruth[it];row.link_class=static_cast<int>(LinkClass::TRUTH_MISS);
    links.push_back(row);
  }
  return links;
}

struct Metadata
{
  std::string schema_sha256,semantic_sha256,source_sha256,model_sha256,config_sha256,code_sha256;
  // Optional for protected historical writers, mandatory in the THE-121/122
  // packet/terminal contract.  This permits old evidence to remain readable
  // while new production binds one campaign-level reconstruction manifest.
  std::string provenance_manifest_sha256;
  std::string analysis_jet_node_suffix;
  std::string data_retention_profile;
  std::string data_retention_contract_sha256;
};

struct CompletionMetadata
{
  std::vector<std::pair<std::string,std::uint64_t>> counters;
};

enum class WriterMode
{
  SERIALIZE,
  VALIDATE_ONLY
};

class Writer
{
 public:
  Writer() = default;
  Writer(const Writer&) = delete;
  Writer& operator=(const Writer&) = delete;

  bool initialize(TFile* file, const Metadata& metadata, std::string* error=nullptr)
  {
    return initialize(file,metadata,WriterMode::SERIALIZE,error);
  }

  bool initialize(TFile* file, const Metadata& metadata, WriterMode mode, std::string* error=nullptr)
  {
    if (!file || !file->IsOpen()) return fail(error,"output file is not open");
    if (!validMetadata(metadata)) return fail(error,"metadata hashes must be 64 lowercase hexadecimal characters");
    if(mode!=WriterMode::SERIALIZE&&mode!=WriterMode::VALIDATE_ONLY)return fail(error,"unsupported writer mode");
    m_file=file; m_metadata=metadata; m_mode=mode;
    m_file->SetCompressionAlgorithm(static_cast<int>(ROOT::RCompressionSetting::EAlgorithm::kZSTD));
    m_file->SetCompressionLevel(5);
    if(m_mode==WriterMode::VALIDATE_ONLY){m_initialized=true;return true;}
    TDirectory* saved=gDirectory;
    m_dir=m_file->GetDirectory("ReplayFoundationV1");
    if (!m_dir) m_dir=m_file->mkdir("ReplayFoundationV1");
    if (!m_dir) return fail(error,"failed to create ReplayFoundationV1 directory");
    m_dir->cd();
    bookTrees();
    TNamed schema("rj_replay_schema",kSchemaName); schema.Write("rj_replay_schema",TObject::kOverwrite);
    const std::string schemaVersionValue=std::to_string(kSchemaVersion);
    TNamed schemaVersion("rj_replay_schema_version",schemaVersionValue.c_str()); schemaVersion.Write("rj_replay_schema_version",TObject::kOverwrite);
    writeMeta("schema_sha256",metadata.schema_sha256); writeMeta("semantic_sha256",metadata.semantic_sha256);
    writeMeta("source_sha256",metadata.source_sha256); writeMeta("model_sha256",metadata.model_sha256);
    writeMeta("config_sha256",metadata.config_sha256); writeMeta("code_sha256",metadata.code_sha256);
    writeMeta("jet_area_definition","fastjet_active_area");
    writeMeta("identity_encoding","SHARD_SCOPED_DOMAIN_AND_ONE_BASED_REFERENCE_V1");
    writeMeta("identity_join_scope","provenance_manifest_sha256,event_id,native_collection_key");
    writeMeta("relationship_validation_scope","EVENT_LOCAL_FOREIGN_KEYS_GLOBAL_EVENT_DUPLICATE_GUARD_V1");
    writeMeta("analysis_jet_node_suffix",metadata.analysis_jet_node_suffix);
    if(!metadata.provenance_manifest_sha256.empty())
      writeMeta("provenance_manifest_sha256",metadata.provenance_manifest_sha256);
    if(!metadata.data_retention_profile.empty())
    {
      writeMeta("data_retention_profile",metadata.data_retention_profile);
      writeMeta("data_retention_contract_sha256",
                metadata.data_retention_contract_sha256);
    }
    if (saved) saved->cd();
    m_initialized=true;
    return true;
  }

  bool fill(const SourceOccurrenceRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!insert(m_sourceIds,r.id))return fail(e,"invalid or duplicate source identity"); const auto sourceRef=ensureRef(m_sourceRefs,r.id); if(validationOnly())return true; m_source=r; serializeRef(m_sourceId,CompactIdentityDomain::SOURCE,sourceRef); m_tSource->Fill(); return true; }
  bool fill(const EventRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_sourceIds,r.source_id)||!validBackgroundWitness(r)||!validSourceOrderWitness(r)||!insert(m_seenEventIds,r.id))return fail(e,"event foreign key, background witness, or identity failure"); resetEventScope(); if(!insert(m_eventIds,r.id))return fail(e,"event scope initialization failure"); const auto eventRef=bindSequentialRef(m_eventRefs,r.id,m_eventRows),sourceRef=lookupRef(m_sourceRefs,r.source_id); if(eventRef==0||sourceRef==0)return fail(e,"event compact reference failure"); if(validationOnly())return true; m_event=r; serializeRef(m_eventId,CompactIdentityDomain::EVENT,eventRef); serializeRef(m_eventSourceId,CompactIdentityDomain::SOURCE,sourceRef); m_eventSequence=serialize(r.event_sequence); m_eventPhysicalEventSequence=serialize(r.physical_event_sequence); m_eventSourceFileOrdinal=serialize(r.source_file_ordinal); m_eventSourceEntryOrdinal=serialize(r.source_entry_ordinal); m_eventSourceGlobalEntryOrdinal=serialize(r.source_global_entry_ordinal); m_eventTriggerBits=serialize(r.trigger_bits); m_eventLiveTriggerBits=serialize(r.live_trigger_bits); m_eventScaledTriggerBits=serialize(r.scaled_trigger_bits); m_eventTriggerPacketStatus=serialize(r.trigger_packet_status); m_eventTriggerScalerSnapshotId=serialize(r.trigger_scaler_snapshot_id); m_eventWeightedFillCount=serialize(r.ppg12_weighted_fill_count); m_eventRawFillCount=serialize(r.ppg12_raw_fill_count); for(std::size_t code=1;code<5;++code){m_eventWeightedFillCountByCode[code]=serialize(r.ppg12_weighted_fill_count_by_code[code]);m_eventRawFillCountByCode[code]=serialize(r.ppg12_raw_fill_count_by_code[code]);} m_tEvent->Fill(); return true; }
  bool fill(const PhotonCandidateRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_candidateIds,r.id))return fail(e,"candidate foreign key or identity failure"); const auto candidateRef=bindSequentialRef(m_candidateRefs,r.id,m_candidateRows),eventRef=lookupRef(m_eventRefs,r.event_id); if(candidateRef==0||eventRef==0)return fail(e,"candidate compact reference failure"); if(validationOnly())return true; m_candidate=r; serializeRef(m_candidateId,CompactIdentityDomain::CANDIDATE,candidateRef); serializeRef(m_candidateEventId,CompactIdentityDomain::EVENT,eventRef); m_candidatePreselectionBitmask=serialize(r.preselection_bitmask); m_tCandidate->Fill(); return true; }
  bool fill(const ModelEvaluationRow& r,std::string* e=nullptr)
  {
    if(!ready(e)||r.model_id.isNull()||r.shower_definition_id.empty()||
       !validHex64(r.shower_semantic_sha256)||!contains(m_candidateIds,r.candidate_id))
      return fail(e,"model evaluation shower semantics or foreign-key failure");
    const auto candidateRef=lookupRef(m_candidateRefs,r.candidate_id);
    const auto modelRef=ensureRef(m_modelRefs,r.model_id);
    if(candidateRef==0||modelRef==0||
       !insert(m_modelEvalIds,Identity128{candidateRef,modelRef}))
      return fail(e,"model compact reference or duplicate failure");
    if(validationOnly())return true;
    m_model=r;
    serializeRef(m_modelCandidateId,CompactIdentityDomain::CANDIDATE,candidateRef);
    serializeRef(m_modelId,CompactIdentityDomain::MODEL,modelRef);
    m_tModel->Fill();
    return true;
  }
  bool fill(const ShowerCellRow& r,std::string* e=nullptr)
  {
    const bool gridOrOwnedProvenance=
        (r.grid_membership_bitmask>0&&r.grid_membership_bitmask<=3)||
        (r.grid_membership_bitmask==0&&r.rawcluster_owned!=0&&
         r.rawcluster_value_present!=0);
    if(!ready(e)||r.tower_eta_index<0||r.tower_eta_index>=96||
       r.tower_phi_index<0||r.tower_phi_index>=256||!gridOrOwnedProvenance||
       !contains(m_candidateIds,r.candidate_id))
      return fail(e,"shower-cell tower identity, provenance, or candidate foreign-key failure");
    const auto candidateRef=lookupRef(m_candidateRefs,r.candidate_id);
    if(candidateRef==0||
       !insert(m_showerCellIds,Identity128{candidateRef,r.tower_key}))
      return fail(e,"shower-cell compact reference or duplicate failure");
    if(validationOnly())return true;
    m_shower=r;
    serializeRef(m_showerCandidateId,CompactIdentityDomain::CANDIDATE,candidateRef);
    m_showerTowerKey=serialize(r.tower_key);
    m_tShower->Fill();
    return true;
  }
  bool fill(const ShowerFeatureViewRow& r,std::string* e=nullptr)
  {
    if(!ready(e)||r.definition_id.isNull()||r.definition_name.empty()||
       !validHex64(r.semantic_sha256)||r.center_eta_index<0||
       r.center_eta_index>=96||r.center_phi_index<0||r.center_phi_index>=256||
       !std::isfinite(r.raw_center_eta)||!std::isfinite(r.raw_center_phi)||
       !contains(m_candidateIds,r.candidate_id))
      return fail(e,"shower-feature-view center, semantic hash, or foreign-key failure");
    const auto candidateRef=lookupRef(m_candidateRefs,r.candidate_id);
    const auto definitionRef=ensureRef(m_showerDefinitionRefs,r.definition_id);
    if(candidateRef==0||definitionRef==0||
       !insert(m_showerViewIds,Identity128{candidateRef,definitionRef}))
      return fail(e,"shower-feature compact reference or duplicate failure");
    if(validationOnly())return true;
    m_showerView=r;
    serializeRef(m_showerViewCandidateId,CompactIdentityDomain::CANDIDATE,candidateRef);
    serializeRef(m_showerViewDefinitionId,CompactIdentityDomain::SHOWER_DEFINITION,definitionRef);
    m_tShowerView->Fill();
    return true;
  }
  bool fill(const IsolationConstituentRow& r,std::string* e=nullptr){ if(!ready(e)||!contains(m_candidateIds,r.candidate_id)||r.native_source<0)return fail(e,"isolation constituent native key or foreign key failure"); const auto candidateRef=lookupRef(m_candidateRefs,r.candidate_id); if(candidateRef==0)return fail(e,"isolation constituent compact reference failure"); const auto rowRef=++m_isoConstituentRows; if(validationOnly())return true; m_isoConstituent=r; serializeRef(m_isoConstituentCandidateId,CompactIdentityDomain::CANDIDATE,candidateRef); serializeRef(m_isoConstituentId,CompactIdentityDomain::ISOLATION_CONSTITUENT,rowRef); m_tIsoConstituent->Fill(); return true; }
  bool fill(const IsolationWitnessRow& r,std::string* e=nullptr){ if(!ready(e)||r.isolation_id.isNull()||!contains(m_candidateIds,r.candidate_id))return fail(e,"isolation witness foreign key failure"); const auto candidateRef=lookupRef(m_candidateRefs,r.candidate_id); if(candidateRef==0)return fail(e,"isolation witness compact reference failure"); const auto rowRef=++m_isoWitnessRows; if(validationOnly())return true; m_isoWitness=r; serializeRef(m_isoWitnessCandidateId,CompactIdentityDomain::CANDIDATE,candidateRef); serializeRef(m_isoWitnessId,CompactIdentityDomain::ISOLATION_WITNESS,rowRef); m_tIsoWitness->Fill(); return true; }
  bool fill(const JetRow& r,std::string* e=nullptr){ const bool areaValid=std::isfinite(r.area)&&r.area>=0.0&&(r.quality_bitmask&kJetFastJetAreaValid)!=0; const bool rawLinkValid=(r.quality_bitmask&kJetCalibRawOrdinalLinked)!=0; if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!areaValid||!rawLinkValid||!insert(m_jetIds,r.id))return fail(e,"jet foreign key, FastJet area, JetCalib raw link, or identity failure"); const auto jetRef=bindSequentialRef(m_jetRefs,r.id,m_jetRows),eventRef=lookupRef(m_eventRefs,r.event_id); if(jetRef==0||eventRef==0)return fail(e,"jet compact reference failure"); if(validationOnly())return true; m_jet=r; serializeRef(m_jetId,CompactIdentityDomain::JET,jetRef); serializeRef(m_jetEventId,CompactIdentityDomain::EVENT,eventRef); m_jetQualityBitmask=serialize(r.quality_bitmask); m_tJet->Fill(); return true; }
  bool fill(const JetConstituentRow& r,std::string* e=nullptr){ if(!ready(e)||!contains(m_jetIds,r.jet_id)||r.native_source<0)return fail(e,"jet constituent native key or foreign key failure"); const auto jetRef=lookupRef(m_jetRefs,r.jet_id); if(jetRef==0)return fail(e,"jet constituent compact reference failure"); const auto rowRef=++m_jetConstituentRows; if(validationOnly())return true; m_jetConstituent=r; serializeRef(m_jetConstituentJetId,CompactIdentityDomain::JET,jetRef); serializeRef(m_jetConstituentId,CompactIdentityDomain::JET_CONSTITUENT,rowRef); m_tJetConstituent->Fill(); return true; }
  bool fill(const PhotonJetPairRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!contains(m_candidateIds,r.candidate_id)||!contains(m_jetIds,r.jet_id)||!insert(m_pairIds,r.id))return fail(e,"pair foreign key or identity failure"); const auto eventRef=lookupRef(m_eventRefs,r.event_id),candidateRef=lookupRef(m_candidateRefs,r.candidate_id),jetRef=lookupRef(m_jetRefs,r.jet_id),rowRef=++m_pairRows; if(eventRef==0||candidateRef==0||jetRef==0)return fail(e,"pair compact reference failure"); if(validationOnly())return true; m_pair=r; serializeRef(m_pairId,CompactIdentityDomain::PAIR,rowRef); serializeRef(m_pairEventId,CompactIdentityDomain::EVENT,eventRef); serializeRef(m_pairCandidateId,CompactIdentityDomain::CANDIDATE,candidateRef); serializeRef(m_pairJetId,CompactIdentityDomain::JET,jetRef); m_tPair->Fill(); return true; }
  bool fill(const TruthPhotonRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_truthPhotonIds,r.id))return fail(e,"truth photon foreign key or identity failure"); const auto truthRef=bindSequentialRef(m_truthPhotonRefs,r.id,m_truthPhotonRows),eventRef=lookupRef(m_eventRefs,r.event_id); if(truthRef==0||eventRef==0)return fail(e,"truth photon compact reference failure"); if(validationOnly())return true; m_truthPhoton=r; serializeRef(m_truthPhotonId,CompactIdentityDomain::TRUTH_PHOTON,truthRef); serializeRef(m_truthPhotonEventId,CompactIdentityDomain::EVENT,eventRef); m_tTruthPhoton->Fill(); return true; }
  bool fill(const TruthPhotonMissOccurrenceRow& r,std::string* e=nullptr){ const bool validPath=r.decision_path=="NO_SELECTED_RECO"||r.decision_path=="SELECTED_RECO_NOT_TRUTH_SIGNAL"; if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||r.trigger.empty()||r.iso_view.empty()||!validPath||!std::isfinite(r.truth_pt)||r.truth_pt<=0.0||r.centrality_bin_index<0||r.occurrence_ordinal<0||!insert(m_truthPhotonMissOccurrenceIds,r.id))return fail(e,"truth photon miss occurrence identity/provenance/foreign-key failure"); const auto eventRef=lookupRef(m_eventRefs,r.event_id),rowRef=++m_occurrenceRows; if(eventRef==0)return fail(e,"truth photon miss compact reference failure"); if(validationOnly())return true; m_truthPhotonMissOccurrence=r; serializeRef(m_truthPhotonMissOccurrenceId,CompactIdentityDomain::OCCURRENCE,rowRef); serializeRef(m_truthPhotonMissOccurrenceEventId,CompactIdentityDomain::EVENT,eventRef); m_tTruthPhotonMissOccurrence->Fill(); return true; }
  bool fill(const EmbeddedPhotonDiagnosticOccurrenceRow& r,std::string* e=nullptr){ const bool audit=r.occurrence_kind=="AUDIT"&&r.audit_bin>=1&&r.audit_bin<=16; const bool photon=r.occurrence_kind=="TRUTH_PHOTON"&&(r.photon_class==1||r.photon_class==2)&&r.centrality_bin_index>=0&&r.centrality_bin_index<=2&&(r.sample_code==12||r.sample_code==20||r.sample_code==30||r.sample_code==40)&&std::isfinite(r.truth_pt)&&r.truth_pt>0.0&&std::isfinite(r.truth_isolation_witness); const bool stitch=r.occurrence_kind=="STITCH"&&r.stitch_decision>=1&&r.stitch_decision<=4&&(r.sample_code==0||r.sample_code==12||r.sample_code==20)&&(!std::isfinite(r.stitch_photon_pt)||r.stitch_photon_pt>0.0); if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||r.trigger.empty()||r.occurrence_ordinal<0||(!audit&&!photon&&!stitch)||!insert(m_embeddedPhotonDiagnosticOccurrenceIds,r.id))return fail(e,"embedded photon diagnostic occurrence identity/provenance/foreign-key failure"); const auto eventRef=lookupRef(m_eventRefs,r.event_id),rowRef=++m_occurrenceRows; if(eventRef==0)return fail(e,"embedded diagnostic compact reference failure"); if(validationOnly())return true; m_embeddedPhotonDiagnosticOccurrence=r; serializeRef(m_embeddedPhotonDiagnosticOccurrenceId,CompactIdentityDomain::OCCURRENCE,rowRef); serializeRef(m_embeddedPhotonDiagnosticOccurrenceEventId,CompactIdentityDomain::EVENT,eventRef); m_tEmbeddedPhotonDiagnosticOccurrence->Fill(); return true; }
  bool fill(const PPG12DiagnosticOccurrenceRow& r,std::string* e=nullptr)
  {
    const bool audit=(r.occurrence_kind=="FIG8_AUDIT"||r.occurrence_kind=="EISO_VERTEX_AUDIT"||r.occurrence_kind=="VERTEX_CONTRACT_AUDIT"||r.occurrence_kind=="FIG3_OWNERSHIP_AUDIT"||r.occurrence_kind=="FIG6_OWNERSHIP_AUDIT")&&r.stage_code>0;
    const bool period=r.occurrence_kind=="PERIOD_CONTRACT_VALUE"&&r.stage_code>=1&&r.stage_code<=13&&std::isfinite(r.flow_value);
    const bool vertex=r.occurrence_kind=="VERTEX_QA_STAGE"&&(r.stage_code==18||r.stage_code==19)&&std::isfinite(r.vertex_weight)&&std::isfinite(r.period_event_weight);
    const bool inclusiveJet=r.occurrence_kind=="PP_INCLUSIVE_JET_EVENT"&&r.sample_code>=1&&r.sample_code<=6&&r.selection_bitmask<=15;
    const bool cluster=r.occurrence_kind=="FIG8_CLUSTER"&&r.cluster_ordinal>=0&&r.truth_track_id>=0&&std::isfinite(r.cluster_energy)&&std::isfinite(r.cluster_eta)&&std::isfinite(r.cluster_et)&&std::isfinite(r.truth_energy)&&std::isfinite(r.truth_eta)&&std::isfinite(r.truth_et)&&std::isfinite(r.energy_contribution);
    const bool stitch=r.occurrence_kind=="PP_STITCH_SOURCE"&&r.source_code>=1&&r.source_code<=7&&r.sample_code>=1&&r.sample_code<=6&&r.decision>=1&&r.decision<=4&&(!std::isfinite(r.source_photon_pt)||r.source_photon_pt>=0.0);
    const bool flow=r.occurrence_kind=="PP_STITCH_FLOW"&&r.stage_code>0&&std::isfinite(r.flow_value);
    const bool fig6=r.occurrence_kind=="PP_FIG6_WEIGHTED_SOURCE"&&r.sample_code>=1&&r.sample_code<=6&&r.decision>=1&&r.decision<=3&&(!std::isfinite(r.flow_value)||r.flow_value>=0.0)&&std::isfinite(r.occurrence_weight);
    const bool photonDecision=(r.occurrence_kind=="PHOTON_ALL_COMMON"||r.occurrence_kind=="PHOTON_TIGHT_ABCD")&&std::isfinite(r.reco_photon_pt)&&r.reco_photon_pt>0.0&&std::isfinite(r.occurrence_weight);
    const bool truthIsolation=r.occurrence_kind=="TRUTH_ISOLATION_FILL"&&
        (r.decision==1||r.decision==2)&&std::isfinite(r.flow_value)&&
        std::isfinite(r.occurrence_weight);
    const bool truthEfficiency=r.occurrence_kind=="TRUTH_EFFICIENCY_DEN_FILL"&&
        std::isfinite(r.truth_photon_pt)&&r.truth_photon_pt>0.0&&
        std::isfinite(r.occurrence_weight);
    const bool fig6Efficiency=r.occurrence_kind=="PP_FIG6_EFFICIENCY_STAGE"&&
        r.stage_code>=1&&r.stage_code<=5&&std::isfinite(r.truth_photon_pt)&&
        r.truth_photon_pt>0.0&&std::isfinite(r.occurrence_weight);
    const bool fig3TruthIsolation=r.occurrence_kind=="PP_FIG3_TRUTH_ISOLATION"&&
        (r.source_code==1||r.source_code==2)&&std::isfinite(r.truth_photon_pt)&&
        r.truth_photon_pt>0.0&&std::isfinite(r.flow_value)&&
        std::isfinite(r.occurrence_weight);
    const bool isolationStack=r.occurrence_kind=="PP_ISOLATION_STACK_FILL"&&
        r.sample_code>=1&&r.sample_code<=11&&(r.decision==1||r.decision==2)&&
        std::isfinite(r.flow_value)&&std::isfinite(r.occurrence_weight);
    const bool truthSpectrum=r.occurrence_kind=="PP_TRUTH_SPECTRUM_FILL"&&
        std::isfinite(r.truth_photon_pt)&&r.truth_photon_pt>0.0&&
        std::isfinite(r.occurrence_weight);
    const bool leadtagTruth=r.occurrence_kind=="PP_LEADTAG_TRUTH_FILL"&&
        (r.decision==1||r.decision==2)&&std::isfinite(r.truth_photon_pt)&&
        r.truth_photon_pt>0.0&&std::isfinite(r.occurrence_weight);
    if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||r.trigger.empty()||
       r.occurrence_ordinal<0||!std::isfinite(r.occurrence_weight)||
       (!audit&&!period&&!vertex&&!inclusiveJet&&!cluster&&!stitch&&!flow&&!fig6&&
        !photonDecision&&!truthIsolation&&!truthEfficiency&&!fig6Efficiency&&
        !fig3TruthIsolation&&!isolationStack&&!truthSpectrum&&!leadtagTruth)||
       !insert(m_ppg12DiagnosticOccurrenceIds,r.id))
      return fail(e,"PPG12 diagnostic occurrence identity/provenance/foreign-key failure");
    const auto eventRef=lookupRef(m_eventRefs,r.event_id),rowRef=++m_occurrenceRows;
    if(eventRef==0)return fail(e,"PPG12 diagnostic compact reference failure");
    if(validationOnly())return true;
    m_ppg12DiagnosticOccurrence=r;
    serializeRef(m_ppg12DiagnosticOccurrenceId,CompactIdentityDomain::OCCURRENCE,rowRef);
    serializeRef(m_ppg12DiagnosticOccurrenceEventId,CompactIdentityDomain::EVENT,eventRef);
    m_ppg12DiagnosticSelectionBitmask=serialize(r.selection_bitmask);
    m_tPPG12DiagnosticOccurrence->Fill();
    return true;
  }
  bool fill(const TruthJetRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_truthJetIds,r.id))return fail(e,"truth jet foreign key or identity failure"); const auto truthRef=bindSequentialRef(m_truthJetRefs,r.id,m_truthJetRows),eventRef=lookupRef(m_eventRefs,r.event_id); if(truthRef==0||eventRef==0)return fail(e,"truth jet compact reference failure"); if(validationOnly())return true; m_truthJet=r; serializeRef(m_truthJetId,CompactIdentityDomain::TRUTH_JET,truthRef); serializeRef(m_truthJetEventId,CompactIdentityDomain::EVENT,eventRef); m_tTruthJet->Fill(); return true; }
  bool fill(const RecoTruthLinkRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!validRecoTruthLink(r)||!insert(m_linkIds,r.id))return fail(e,"reco-truth link identity/type/foreign-key failure"); const auto rowRef=++m_linkRows; if(validationOnly())return true; m_link=r; serializeRef(m_linkId,CompactIdentityDomain::LINK,rowRef); if(!serializeTypedRecoRef(m_linkRecoId,r.reco_type,r.reco_id)||!serializeTypedTruthRef(m_linkTruthId,r.truth_type,r.truth_id))return fail(e,"reco-truth compact reference failure"); m_tLink->Fill(); return true; }
  bool fill(const WeightComponentRow& r,std::string* e=nullptr){ if(!ready(e)||r.target_id.isNull()||r.component_type.empty()||r.application_count<0||!validWeightTarget(r))return fail(e,"weight reference/type/count failure"); if(validationOnly())return true; m_weight=r; if(!serializeWeightTarget(m_weightTargetId,r))return fail(e,"weight compact reference failure"); m_tWeight->Fill(); return true; }
  bool fill(const EventDisplaySnapshotRow& r,std::string* e=nullptr){ if(!ready(e)||r.id.isNull()||!contains(m_eventIds,r.event_id)||!insert(m_snapshotIds,r.id))return fail(e,"snapshot foreign key or identity failure"); const auto eventRef=lookupRef(m_eventRefs,r.event_id),rowRef=++m_snapshotRows; if(eventRef==0)return fail(e,"snapshot compact reference failure"); if(validationOnly())return true; m_snapshot=r; serializeRef(m_snapshotId,CompactIdentityDomain::SNAPSHOT,rowRef); serializeRef(m_snapshotEventId,CompactIdentityDomain::EVENT,eventRef); m_tSnapshot->Fill(); return true; }

  bool finish(std::string* error=nullptr)
  {
    return finish(CompletionMetadata{},error);
  }

  bool finish(const CompletionMetadata& completion,std::string* error=nullptr)
  {
    if(!ready(error)) return false;
    if(validationOnly()){m_finished=true;return true;}
    TDirectory* saved=gDirectory; m_dir->cd();
    for(TTree* tree:m_trees) if(!tree||tree->Write("",TObject::kOverwrite)<=0){ if(saved)saved->cd(); return fail(error,"tree write failure"); }
    for(const auto& item:completion.counters)
      writeMeta(item.first.c_str(),std::to_string(item.second));
    TNamed complete("rj_replay_complete","1"); complete.Write("rj_replay_complete",TObject::kOverwrite);
    if(saved) saved->cd();
    m_finished=true;
    return true;
  }

  const Metadata& metadata() const { return m_metadata; }

 private:
  static bool fail(std::string* error,const std::string& message){ if(error)*error=message; return false; }
  bool ready(std::string* e) const { return m_initialized&&!m_finished?true:fail(e,"writer is not active"); }
  bool validationOnly() const { return m_mode==WriterMode::VALIDATE_ONLY; }
  static bool validHex64(const std::string& s){ if(s.size()!=64)return false; for(char c:s)if(!((c>='0'&&c<='9')||(c>='a'&&c<='f')))return false; return true; }
  static bool validMetadata(const Metadata& m){ return validHex64(m.schema_sha256)&&validHex64(m.semantic_sha256)&&validHex64(m.source_sha256)&&validHex64(m.model_sha256)&&validHex64(m.config_sha256)&&validHex64(m.code_sha256)&&(m.provenance_manifest_sha256.empty()||validHex64(m.provenance_manifest_sha256))&&validAnalysisJetNodeSuffix(m.analysis_jet_node_suffix)&&(m.data_retention_profile.empty()||validHex64(m.data_retention_contract_sha256)); }
  static bool finiteVector(const std::vector<float>& values){ return std::all_of(values.begin(),values.end(),[](float value){return std::isfinite(value);}); }
  static bool validBackgroundWitness(const EventRow& r)
  {
    if(r.jet_background_contract_state==0)
      return r.jet_background_flow_mode==-1&&r.jet_background_eta_bins==0&&
             !std::isfinite(r.jet_background_v2)&&!std::isfinite(r.jet_background_psi2)&&
             r.jet_background_ue_emcal.empty()&&r.jet_background_ue_ihcal.empty()&&
             r.jet_background_ue_ohcal.empty()&&r.jet_background_n_strips==-1&&
             r.jet_background_n_towers==-1&&r.jet_background_flow_failure_state==-1;
    const std::size_t bins=r.jet_background_ue_emcal.size();
    return r.jet_background_contract_state==1&&
           r.jet_background_flow_mode>=0&&r.jet_background_flow_mode<=3&&bins>0&&
           r.jet_background_ue_ihcal.size()==bins&&r.jet_background_ue_ohcal.size()==bins&&
           r.jet_background_eta_bins==static_cast<std::int32_t>(bins)&&
           std::isfinite(r.jet_background_v2)&&std::isfinite(r.jet_background_psi2)&&
           finiteVector(r.jet_background_ue_emcal)&&finiteVector(r.jet_background_ue_ihcal)&&
           finiteVector(r.jet_background_ue_ohcal)&&r.jet_background_n_strips>=0&&
           r.jet_background_n_towers>=0&&
           (r.jet_background_flow_failure_state==0||r.jet_background_flow_failure_state==1);
  }
  static bool validSourceOrderWitness(const EventRow& r)
  {
    const bool absent=r.source_file_ordinal==-1&&r.source_entry_ordinal==-1&&
        r.source_global_entry_ordinal==-1;
    const bool present=r.source_file_ordinal>=0&&r.source_entry_ordinal>=0&&
        r.source_global_entry_ordinal>=r.source_entry_ordinal;
    return absent||present;
  }
  static bool insert(std::unordered_set<Identity128,IdentityHash>& s,const Identity128& id){ return s.insert(id).second; }
  static bool contains(const std::unordered_set<Identity128,IdentityHash>& s,const Identity128& id){ return s.find(id)!=s.end(); }
  using IdentityRefMap=std::unordered_map<Identity128,std::uint64_t,IdentityHash>;
  static std::uint64_t ensureRef(IdentityRefMap& refs,const Identity128& id)
  {
    if(id.isNull())return 0;
    const auto found=refs.find(id);
    if(found!=refs.end())return found->second;
    const std::uint64_t ref=static_cast<std::uint64_t>(refs.size())+1ULL;
    refs.emplace(id,ref);
    return ref;
  }
  static std::uint64_t lookupRef(const IdentityRefMap& refs,const Identity128& id)
  {
    const auto found=refs.find(id);
    return found==refs.end()?0ULL:found->second;
  }
  static std::uint64_t bindSequentialRef(IdentityRefMap& refs,const Identity128& id,std::uint64_t& counter)
  {
    if(id.isNull()||refs.find(id)!=refs.end())return 0;
    const std::uint64_t ref=++counter;
    refs.emplace(id,ref);
    return ref;
  }
  void resetEventScope()
  {
    m_eventIds.clear(); m_candidateIds.clear(); m_modelEvalIds.clear();
    m_showerCellIds.clear(); m_showerViewIds.clear(); m_jetIds.clear();
    m_pairIds.clear(); m_truthPhotonIds.clear();
    m_truthPhotonMissOccurrenceIds.clear();
    m_embeddedPhotonDiagnosticOccurrenceIds.clear();
    m_ppg12DiagnosticOccurrenceIds.clear(); m_truthJetIds.clear();
    m_linkIds.clear(); m_snapshotIds.clear();
    m_eventRefs.clear(); m_candidateRefs.clear(); m_jetRefs.clear();
    m_truthPhotonRefs.clear(); m_truthJetRefs.clear();
  }
  static SerializedUInt64 serialize(std::uint64_t value){ return static_cast<SerializedUInt64>(value); }
  static SerializedInt64 serialize(std::int64_t value){ return static_cast<SerializedInt64>(value); }
  static void serialize(SerializedIdentity128& output,const Identity128& input)
  { output.hi=serialize(input.hi); output.lo=serialize(input.lo); }
  static void serializeRef(SerializedIdentity128& output,CompactIdentityDomain domain,std::uint64_t ref)
  {
    output.hi=serialize(static_cast<std::uint64_t>(domain));
    output.lo=serialize(ref);
  }
  static void serializeNullRef(SerializedIdentity128& output){output.hi=0;output.lo=0;}
  bool serializeTypedRecoRef(SerializedIdentity128& output,std::int32_t type,const Identity128& id) const
  {
    const auto typed=static_cast<RecoTruthType>(type);
    if(typed==RecoTruthType::NONE){serializeNullRef(output);return id.isNull();}
    if(typed==RecoTruthType::PHOTON){const auto ref=lookupRef(m_candidateRefs,id);if(ref==0)return false;serializeRef(output,CompactIdentityDomain::CANDIDATE,ref);return true;}
    if(typed==RecoTruthType::JET){const auto ref=lookupRef(m_jetRefs,id);if(ref==0)return false;serializeRef(output,CompactIdentityDomain::JET,ref);return true;}
    return false;
  }
  bool serializeTypedTruthRef(SerializedIdentity128& output,std::int32_t type,const Identity128& id) const
  {
    const auto typed=static_cast<RecoTruthType>(type);
    if(typed==RecoTruthType::NONE){serializeNullRef(output);return id.isNull();}
    if(typed==RecoTruthType::PHOTON){const auto ref=lookupRef(m_truthPhotonRefs,id);if(ref==0)return false;serializeRef(output,CompactIdentityDomain::TRUTH_PHOTON,ref);return true;}
    if(typed==RecoTruthType::JET){const auto ref=lookupRef(m_truthJetRefs,id);if(ref==0)return false;serializeRef(output,CompactIdentityDomain::TRUTH_JET,ref);return true;}
    return false;
  }
  bool validWeightTarget(const WeightComponentRow& row) const
  {
    switch(static_cast<WeightTargetType>(row.target_type))
    {
      case WeightTargetType::EVENT:return contains(m_eventIds,row.target_id);
      case WeightTargetType::CANDIDATE:return contains(m_candidateIds,row.target_id);
      case WeightTargetType::JET:return contains(m_jetIds,row.target_id);
      case WeightTargetType::TRUTH_PHOTON:return contains(m_truthPhotonIds,row.target_id);
      case WeightTargetType::TRUTH_JET:return contains(m_truthJetIds,row.target_id);
    }
    return false;
  }
  bool serializeWeightTarget(SerializedIdentity128& output,const WeightComponentRow& row) const
  {
    std::uint64_t ref=0;
    CompactIdentityDomain domain=CompactIdentityDomain::NULL_REF;
    switch(static_cast<WeightTargetType>(row.target_type))
    {
      case WeightTargetType::EVENT:ref=lookupRef(m_eventRefs,row.target_id);domain=CompactIdentityDomain::EVENT;break;
      case WeightTargetType::CANDIDATE:ref=lookupRef(m_candidateRefs,row.target_id);domain=CompactIdentityDomain::CANDIDATE;break;
      case WeightTargetType::JET:ref=lookupRef(m_jetRefs,row.target_id);domain=CompactIdentityDomain::JET;break;
      case WeightTargetType::TRUTH_PHOTON:ref=lookupRef(m_truthPhotonRefs,row.target_id);domain=CompactIdentityDomain::TRUTH_PHOTON;break;
      case WeightTargetType::TRUTH_JET:ref=lookupRef(m_truthJetRefs,row.target_id);domain=CompactIdentityDomain::TRUTH_JET;break;
    }
    if(ref==0)return false;
    serializeRef(output,domain,ref);
    return true;
  }
  void writeMeta(const char* key,const std::string& value){ TNamed obj(key,value.c_str()); obj.Write(key,TObject::kOverwrite); }
  TTree* tree(const char* name){ auto* t=new TTree(name,name); t->SetAutoFlush(-5000000); m_trees.push_back(t); return t; }
  static void unsigned64(TTree* t,const char* n,SerializedUInt64* p){ const std::string leaf=std::string(n)+"/l"; t->Branch(n,p,leaf.c_str()); }
  static void signed64(TTree* t,const char* n,SerializedInt64* p){ const std::string leaf=std::string(n)+"/L"; t->Branch(n,p,leaf.c_str()); }
  static void idBranches(TTree* t,const char* prefix,SerializedIdentity128* id){ const std::string hi=std::string(prefix)+"_hi"; const std::string lo=std::string(prefix)+"_lo"; unsigned64(t,hi.c_str(),&id->hi); unsigned64(t,lo.c_str(),&id->lo); }
  template<class T> static void scalar(TTree* t,const char* n,T* p,const char* type){ const std::string leaf=std::string(n)+"/"+type; t->Branch(n,p,leaf.c_str()); }
  static void str(TTree* t,const char* n,std::string* p){ t->Branch(n,p); }
  template<class T> static void vec(TTree* t,const char* n,std::vector<T>* p){ t->Branch(n,p); }

  void bookTrees()
  {
    m_tSource=tree("RJSourceOccurrenceV1"); idBranches(m_tSource,"source_occurrence_id",&m_sourceId); str(m_tSource,"lane",&m_source.lane); str(m_tSource,"dataset",&m_source.dataset); str(m_tSource,"sample",&m_source.sample); str(m_tSource,"period",&m_source.period); scalar(m_tSource,"run",&m_source.run,"I"); scalar(m_tSource,"segment",&m_source.segment,"I"); str(m_tSource,"si_di_role",&m_source.si_di_role); str(m_tSource,"ownership_state",&m_source.ownership_state); str(m_tSource,"input_uri_hash",&m_source.input_uri_hash); str(m_tSource,"input_file_sha256",&m_source.input_file_sha256); str(m_tSource,"source_manifest_sha256",&m_source.source_manifest_sha256);
    m_tEvent=tree("RJEventV1"); idBranches(m_tEvent,"event_id",&m_eventId); idBranches(m_tEvent,"source_occurrence_id",&m_eventSourceId); scalar(m_tEvent,"run",&m_event.run,"I"); signed64(m_tEvent,"event_sequence",&m_eventSequence); signed64(m_tEvent,"physical_event_sequence",&m_eventPhysicalEventSequence); unsigned64(m_tEvent,"trigger_bits",&m_eventTriggerBits); unsigned64(m_tEvent,"live_trigger_bits",&m_eventLiveTriggerBits); unsigned64(m_tEvent,"scaled_trigger_bits",&m_eventScaledTriggerBits); unsigned64(m_tEvent,"trigger_packet_status",&m_eventTriggerPacketStatus); scalar(m_tEvent,"trigger_capture_version",&m_event.trigger_capture_version,"I"); scalar(m_tEvent,"trigger_decision_state",&m_event.trigger_decision_state,"I"); scalar(m_tEvent,"trigger_decision_available",&m_event.trigger_decision_available,"i"); scalar(m_tEvent,"trigger_packet_version",&m_event.trigger_packet_version,"I"); unsigned64(m_tEvent,"trigger_scaler_snapshot_id",&m_eventTriggerScalerSnapshotId); scalar(m_tEvent,"vertex_z",&m_event.vertex_z,"D"); scalar(m_tEvent,"centrality",&m_event.centrality,"D"); scalar(m_tEvent,"event_weight",&m_event.event_weight,"D"); scalar(m_tEvent,"truth_vertex_z",&m_event.truth_vertex_z,"D"); scalar(m_tEvent,"truth_mb_vertex_z",&m_event.truth_mb_vertex_z,"D"); scalar(m_tEvent,"mbd_total_charge",&m_event.mbd_total_charge,"D"); scalar(m_tEvent,"centrality_replay_version",&m_event.centrality_replay_version,"I"); scalar(m_tEvent,"centrality_pmt_available",&m_event.centrality_pmt_available,"I"); scalar(m_tEvent,"centrality_inputs_valid",&m_event.centrality_inputs_valid,"I"); scalar(m_tEvent,"centrality_mb_decision",&m_event.centrality_mb_decision,"I"); scalar(m_tEvent,"centrality_mbd_z_valid",&m_event.centrality_mbd_z_valid,"I"); scalar(m_tEvent,"centrality_mbd_event",&m_event.centrality_mbd_event,"I"); scalar(m_tEvent,"centrality_mbd_clock",&m_event.centrality_mbd_clock,"I"); scalar(m_tEvent,"centrality_mbd_femclock",&m_event.centrality_mbd_femclock,"I"); scalar(m_tEvent,"centrality_native_bin",&m_event.centrality_native_bin,"I"); scalar(m_tEvent,"centrality_native_valid",&m_event.centrality_native_valid,"I"); scalar(m_tEvent,"centrality_mbd_z",&m_event.centrality_mbd_z,"D"); scalar(m_tEvent,"centrality_selected_charge",&m_event.centrality_selected_charge,"D"); scalar(m_tEvent,"centrality_native_centile",&m_event.centrality_native_centile,"D"); vec(m_tEvent,"centrality_pmt_id",&m_event.centrality_pmt_id); vec(m_tEvent,"centrality_pmt_charge",&m_event.centrality_pmt_charge); vec(m_tEvent,"centrality_pmt_time",&m_event.centrality_pmt_time); vec(m_tEvent,"centrality_pmt_valid",&m_event.centrality_pmt_valid); vec(m_tEvent,"centrality_pmt_selected",&m_event.centrality_pmt_selected); scalar(m_tEvent,"collaborator_interface_version",&m_event.collaborator_interface_version,"I"); scalar(m_tEvent,"reco_vertex_valid",&m_event.reco_vertex_valid,"I"); scalar(m_tEvent,"sample_weight_valid",&m_event.sample_weight_valid,"I"); scalar(m_tEvent,"vertex_weight_valid",&m_event.vertex_weight_valid,"I"); scalar(m_tEvent,"truth_denominator_complete",&m_event.truth_denominator_complete,"I"); scalar(m_tEvent,"mbd_pmt_available",&m_event.mbd_pmt_available,"I"); scalar(m_tEvent,"legacy_event_weight",&m_event.legacy_event_weight,"D"); vec(m_tEvent,"mbd_pmt_id",&m_event.mbd_pmt_id); vec(m_tEvent,"mbd_pmt_arm",&m_event.mbd_pmt_arm); vec(m_tEvent,"mbd_pmt_charge",&m_event.mbd_pmt_charge); vec(m_tEvent,"mbd_pmt_charge_valid",&m_event.mbd_pmt_charge_valid); scalar(m_tEvent,"emcal_total_energy",&m_event.emcal_total_energy,"D"); scalar(m_tEvent,"ihcal_total_energy",&m_event.ihcal_total_energy,"D"); scalar(m_tEvent,"ohcal_total_energy",&m_event.ohcal_total_energy,"D"); scalar(m_tEvent,"total_calo_energy",&m_event.total_calo_energy,"D"); scalar(m_tEvent,"max_cluster_energy",&m_event.max_cluster_energy,"D"); scalar(m_tEvent,"event_calo_require_isgood",&m_event.event_calo_require_isgood,"I"); scalar(m_tEvent,"embedded_minbias_decision",&m_event.embedded_minbias_decision,"I"); scalar(m_tEvent,"jet_background_contract_state",&m_event.jet_background_contract_state,"I"); scalar(m_tEvent,"jet_background_flow_mode",&m_event.jet_background_flow_mode,"I"); scalar(m_tEvent,"jet_background_eta_bins",&m_event.jet_background_eta_bins,"I"); scalar(m_tEvent,"jet_background_v2",&m_event.jet_background_v2,"D"); scalar(m_tEvent,"jet_background_psi2",&m_event.jet_background_psi2,"D"); vec(m_tEvent,"jet_background_ue_emcal",&m_event.jet_background_ue_emcal); vec(m_tEvent,"jet_background_ue_ihcal",&m_event.jet_background_ue_ihcal); vec(m_tEvent,"jet_background_ue_ohcal",&m_event.jet_background_ue_ohcal); scalar(m_tEvent,"jet_background_n_strips",&m_event.jet_background_n_strips,"I"); scalar(m_tEvent,"jet_background_n_towers",&m_event.jet_background_n_towers,"I"); scalar(m_tEvent,"jet_background_flow_failure_state",&m_event.jet_background_flow_failure_state,"I"); scalar(m_tEvent,"nodes_ready",&m_event.nodes_ready,"I"); scalar(m_tEvent,"ppg12_weight_active_state",&m_event.ppg12_weight_active_state,"I"); scalar(m_tEvent,"ppg12_weight_lane_component_code",&m_event.ppg12_weight_lane_component_code,"I"); signed64(m_tEvent,"ppg12_weighted_fill_count",&m_eventWeightedFillCount); signed64(m_tEvent,"ppg12_raw_fill_count",&m_eventRawFillCount); for(std::size_t code=1;code<5;++code){const std::string suffix=std::to_string(code);const std::string weightedName="ppg12_weighted_fill_count_code"+suffix;const std::string rawName="ppg12_raw_fill_count_code"+suffix;signed64(m_tEvent,weightedName.c_str(),&m_eventWeightedFillCountByCode[code]);signed64(m_tEvent,rawName.c_str(),&m_eventRawFillCountByCode[code]);} scalar(m_tEvent,"terminal_status",&m_event.terminal_status,"I"); scalar(m_tEvent,"candidate_count",&m_event.candidate_count,"I"); scalar(m_tEvent,"tag_count",&m_event.tag_count,"I"); scalar(m_tEvent,"recoil_count",&m_event.recoil_count,"I");
    signed64(m_tEvent,"source_file_ordinal",&m_eventSourceFileOrdinal);
    signed64(m_tEvent,"source_entry_ordinal",&m_eventSourceEntryOrdinal);
    signed64(m_tEvent,"source_global_entry_ordinal",&m_eventSourceGlobalEntryOrdinal);
    scalar(m_tEvent,"ppg12_source_arbitration_state",&m_event.ppg12_source_arbitration_state,"I");
    scalar(m_tEvent,"ppg12_leading_source_cluster_key",&m_event.ppg12_leading_source_cluster_key,"i");
    scalar(m_tEvent,"ppg12_leading_source_encounter_ordinal",&m_event.ppg12_leading_source_encounter_ordinal,"I");
    scalar(m_tEvent,"ppg12_leading_source_et",&m_event.ppg12_leading_source_et,"D");
    scalar(m_tEvent,"ppg12_leading_source_eta",&m_event.ppg12_leading_source_eta,"D");
    scalar(m_tEvent,"ppg12_leading_source_phi",&m_event.ppg12_leading_source_phi,"D");
    m_tCandidate=tree("RJPhotonCandidateV1"); idBranches(m_tCandidate,"candidate_id",&m_candidateId); idBranches(m_tCandidate,"event_id",&m_candidateEventId); scalar(m_tCandidate,"encounter_ordinal",&m_candidate.encounter_ordinal,"I"); vec(m_tCandidate,"rank_keys",&m_candidate.rank_keys); scalar(m_tCandidate,"cluster_et",&m_candidate.cluster_et,"D"); scalar(m_tCandidate,"eta",&m_candidate.eta,"D"); scalar(m_tCandidate,"phi",&m_candidate.phi,"D"); scalar(m_tCandidate,"native_weta_cogx",&m_candidate.native_weta_cogx,"D"); scalar(m_tCandidate,"native_wphi_cogx",&m_candidate.native_wphi_cogx,"D"); scalar(m_tCandidate,"native_weta33_cogx",&m_candidate.native_weta33_cogx,"D"); scalar(m_tCandidate,"native_wphi33_cogx",&m_candidate.native_wphi33_cogx,"D"); scalar(m_tCandidate,"native_weta35_cogx",&m_candidate.native_weta35_cogx,"D"); scalar(m_tCandidate,"native_wphi53_cogx",&m_candidate.native_wphi53_cogx,"D"); scalar(m_tCandidate,"native_et1",&m_candidate.native_et1,"D"); scalar(m_tCandidate,"native_e11_over_e33",&m_candidate.native_e11_over_e33,"D"); scalar(m_tCandidate,"native_e32_over_e35",&m_candidate.native_e32_over_e35,"D"); scalar(m_tCandidate,"npb_score",&m_candidate.npb_score,"D"); scalar(m_tCandidate,"auau_npb_score",&m_candidate.auau_npb_score,"D"); vec(m_tCandidate,"ordered_features",&m_candidate.ordered_features); vec(m_tCandidate,"shower_definition_views",&m_candidate.shower_definition_views); unsigned64(m_tCandidate,"preselection_bitmask",&m_candidatePreselectionBitmask); scalar(m_tCandidate,"finite_feature_state",&m_candidate.finite_feature_state,"I"); scalar(m_tCandidate,"below15_retention_state",&m_candidate.below15_retention_state,"I"); scalar(m_tCandidate,"reference_preselection_state",&m_candidate.reference_preselection_state,"I"); scalar(m_tCandidate,"active_preselection_state",&m_candidate.active_preselection_state,"I"); scalar(m_tCandidate,"truth_signal_match_state",&m_candidate.truth_signal_match_state,"I"); scalar(m_tCandidate,"truth_signal_match_barcode",&m_candidate.truth_signal_match_barcode,"I");
    scalar(m_tCandidate,"ppg12_tower_mask_state",&m_candidate.ppg12_tower_mask_state,"I");
    scalar(m_tCandidate,"ppg12_source_eligible_state",&m_candidate.ppg12_source_eligible_state,"I");
    scalar(m_tEvent,"physical_event_sequence_valid",&m_event.physical_event_sequence_valid,"I");
    scalar(m_tEvent,"truth_photon_capture_state",&m_event.truth_photon_inventory.capture_state,"I");
    scalar(m_tEvent,"truth_photon_native_embedded_primary_pid22_count",&m_event.truth_photon_inventory.native_embedded_primary_pid22_count,"I");
    scalar(m_tEvent,"truth_photon_serialized_raw_count",&m_event.truth_photon_inventory.serialized_raw_count,"I");
    scalar(m_tEvent,"truth_photon_rejected_count",&m_event.truth_photon_inventory.rejected_count,"I");
    scalar(m_tEvent,"truth_photon_rejected_invalid_kinematics_count",&m_event.truth_photon_inventory.rejected_invalid_kinematics_count,"I");
    scalar(m_tEvent,"truth_photon_rejected_nonfinite_isolation_count",&m_event.truth_photon_inventory.rejected_nonfinite_isolation_count,"I");
    scalar(m_tEvent,"truth_photon_duplicate_track_count",&m_event.truth_photon_inventory.duplicate_track_count,"I");
    scalar(m_tEvent,"truth_photon_isolation_input_incomplete_count",&m_event.truth_photon_inventory.isolation_input_incomplete_count,"I");
    scalar(m_tEvent,"truth_photon_analysis_signal_count",&m_event.truth_photon_inventory.analysis_signal_count,"I");
    scalar(m_tEvent,"truth_vertex_capture_state",&m_event.truth_vertex_capture_state,"I");
    scalar(m_tEvent,"truth_primary_vertex_id",&m_event.truth_primary_vertex_id,"I");
    vec(m_tEvent,"truth_vertex_id",&m_event.truth_vertex_id);
    vec(m_tEvent,"truth_vertex_embedding_id",&m_event.truth_vertex_embedding_id);
    vec(m_tEvent,"truth_vertex_embedding_valid",&m_event.truth_vertex_embedding_valid);
    vec(m_tEvent,"truth_vertex_z_values",&m_event.truth_vertex_z_values);
    vec(m_tEvent,"truth_vertex_z_valid",&m_event.truth_vertex_z_valid);
    scalar(m_tEvent,"truth_jet_capture_state",&m_event.truth_jet_capture_state,"I");
    vec(m_tEvent,"truth_jet_radius_code",&m_event.truth_jet_radius_code);
    vec(m_tEvent,"truth_jet_container_valid",&m_event.truth_jet_container_valid);
    vec(m_tEvent,"reco_jet_view_radius_code",&m_event.reco_jet_view_radius_code);
    vec(m_tEvent,"reco_jet_view_available",&m_event.reco_jet_view_available);
    vec(m_tEvent,"reco_jet_view_input_identity",&m_event.reco_jet_view_input_identity);
    vec(m_tEvent,"reco_jet_view_subtraction_identity",&m_event.reco_jet_view_subtraction_identity);
    scalar(m_tCandidate,"native_cluster_key",&m_candidate.native_cluster_key,"i");
    scalar(m_tCandidate,"dominant_truth_state",&m_candidate.dominant_truth.state,"I");
    scalar(m_tCandidate,"dominant_truth_evaluator_mode",&m_candidate.dominant_truth.evaluator_mode,"I");
    scalar(m_tCandidate,"dominant_truth_track_id",&m_candidate.dominant_truth.track_id,"I");
    scalar(m_tCandidate,"dominant_truth_pid",&m_candidate.dominant_truth.pid,"I");
    scalar(m_tCandidate,"dominant_truth_barcode",&m_candidate.dominant_truth.barcode,"I");
    scalar(m_tCandidate,"dominant_truth_embedding_id",&m_candidate.dominant_truth.embedding_id,"I");
    scalar(m_tCandidate,"dominant_truth_energy_contribution",&m_candidate.dominant_truth.energy_contribution,"D");
    scalar(m_tCandidate,"dominant_truth_vertex_id",&m_candidate.dominant_truth.vertex_id,"I");
    m_tModel=tree("RJModelEvaluationV1"); idBranches(m_tModel,"candidate_id",&m_modelCandidateId); idBranches(m_tModel,"model_id",&m_modelId); str(m_tModel,"model_sha256",&m_model.model_sha256); str(m_tModel,"shower_definition_id",&m_model.shower_definition_id); str(m_tModel,"shower_semantic_sha256",&m_model.shower_semantic_sha256); vec(m_tModel,"ordered_input_witnesses",&m_model.ordered_input_witnesses); scalar(m_tModel,"raw_score",&m_model.raw_score,"D"); scalar(m_tModel,"finite_score",&m_model.finite_score,"I"); scalar(m_tModel,"applicability_state",&m_model.applicability_state,"I"); scalar(m_tModel,"wp70",&m_model.wp70,"D"); scalar(m_tModel,"wp80",&m_model.wp80,"D"); scalar(m_tModel,"wp90",&m_model.wp90,"D"); scalar(m_tModel,"delta_wp70",&m_model.delta_wp70,"D"); scalar(m_tModel,"delta_wp80",&m_model.delta_wp80,"D"); scalar(m_tModel,"delta_wp90",&m_model.delta_wp90,"D");
    m_tShower=tree("RJShowerCellV1"); idBranches(m_tShower,"candidate_id",&m_showerCandidateId); scalar(m_tShower,"local_eta_index",&m_shower.local_eta_index,"I"); scalar(m_tShower,"local_phi_index",&m_shower.local_phi_index,"I"); scalar(m_tShower,"tower_eta_index",&m_shower.tower_eta_index,"I"); scalar(m_tShower,"tower_phi_index",&m_shower.tower_phi_index,"I"); unsigned64(m_tShower,"tower_key",&m_showerTowerKey); scalar(m_tShower,"calibrated_energy",&m_shower.calibrated_energy,"D"); scalar(m_tShower,"rawcluster_map_value",&m_shower.rawcluster_map_value,"D"); scalar(m_tShower,"is_good",&m_shower.is_good,"I"); scalar(m_tShower,"is_zero",&m_shower.is_zero,"I"); scalar(m_tShower,"is_negative",&m_shower.is_negative,"I"); scalar(m_tShower,"is_nonfinite",&m_shower.is_nonfinite,"I"); scalar(m_tShower,"seed_state",&m_shower.seed_state,"I"); scalar(m_tShower,"denominator_membership",&m_shower.denominator_membership,"I"); scalar(m_tShower,"rawcluster_owned",&m_shower.rawcluster_owned,"I"); scalar(m_tShower,"rawcluster_value_present",&m_shower.rawcluster_value_present,"I"); scalar(m_tShower,"floor0_membership",&m_shower.floor0_membership,"I"); scalar(m_tShower,"floor70_membership",&m_shower.floor70_membership,"I"); scalar(m_tShower,"grid_membership_bitmask",&m_shower.grid_membership_bitmask,"I");
    m_tShowerView=tree("RJShowerFeatureViewV1"); idBranches(m_tShowerView,"candidate_id",&m_showerViewCandidateId); idBranches(m_tShowerView,"definition_id",&m_showerViewDefinitionId); str(m_tShowerView,"definition_name",&m_showerView.definition_name); str(m_tShowerView,"semantic_sha256",&m_showerView.semantic_sha256); vec(m_tShowerView,"ordered_features",&m_showerView.ordered_features); scalar(m_tShowerView,"floor_gev",&m_showerView.floor_gev,"D"); scalar(m_tShowerView,"cog_eta",&m_showerView.cog_eta,"D"); scalar(m_tShowerView,"cog_phi",&m_showerView.cog_phi,"D"); scalar(m_tShowerView,"raw_center_eta",&m_showerView.raw_center_eta,"D"); scalar(m_tShowerView,"raw_center_phi",&m_showerView.raw_center_phi,"D"); scalar(m_tShowerView,"center_eta_index",&m_showerView.center_eta_index,"I"); scalar(m_tShowerView,"center_phi_index",&m_showerView.center_phi_index,"I"); scalar(m_tShowerView,"e11",&m_showerView.e11,"D"); scalar(m_tShowerView,"e33",&m_showerView.e33,"D"); scalar(m_tShowerView,"e32",&m_showerView.e32,"D"); scalar(m_tShowerView,"e35",&m_showerView.e35,"D"); scalar(m_tShowerView,"e11_over_e33",&m_showerView.e11_over_e33,"D"); scalar(m_tShowerView,"e32_over_e35",&m_showerView.e32_over_e35,"D"); scalar(m_tShowerView,"weta_cogx",&m_showerView.weta_cogx,"D"); scalar(m_tShowerView,"wphi_cogx",&m_showerView.wphi_cogx,"D"); scalar(m_tShowerView,"weta33_cogx",&m_showerView.weta33_cogx,"D"); scalar(m_tShowerView,"wphi33_cogx",&m_showerView.wphi33_cogx,"D"); scalar(m_tShowerView,"native_et1",&m_showerView.native_et1,"D"); scalar(m_tShowerView,"native_et2",&m_showerView.native_et2,"D"); scalar(m_tShowerView,"native_et3",&m_showerView.native_et3,"D"); scalar(m_tShowerView,"native_et4",&m_showerView.native_et4,"D"); scalar(m_tShowerView,"moment_eta_numerator",&m_showerView.moment_eta_numerator,"D"); scalar(m_tShowerView,"moment_phi_numerator",&m_showerView.moment_phi_numerator,"D"); scalar(m_tShowerView,"moment_denominator",&m_showerView.moment_denominator,"D"); scalar(m_tShowerView,"moment33_eta_numerator",&m_showerView.moment33_eta_numerator,"D"); scalar(m_tShowerView,"moment33_phi_numerator",&m_showerView.moment33_phi_numerator,"D"); scalar(m_tShowerView,"moment33_denominator",&m_showerView.moment33_denominator,"D"); scalar(m_tShowerView,"energy_source",&m_showerView.energy_source,"I"); scalar(m_tShowerView,"rectangular_membership",&m_showerView.rectangular_membership,"I"); scalar(m_tShowerView,"moment_membership",&m_showerView.moment_membership,"I"); scalar(m_tShowerView,"finite_feature_state",&m_showerView.finite_feature_state,"I"); scalar(m_tShowerView,"good_cell_count",&m_showerView.good_cell_count,"I"); scalar(m_tShowerView,"owned_cell_count",&m_showerView.owned_cell_count,"I"); scalar(m_tShowerView,"active_sum_cell_count",&m_showerView.active_sum_cell_count,"I"); scalar(m_tShowerView,"active_moment_cell_count",&m_showerView.active_moment_cell_count,"I"); scalar(m_tShowerView,"exact_zero_count",&m_showerView.exact_zero_count,"I"); scalar(m_tShowerView,"negative_count",&m_showerView.negative_count,"I"); scalar(m_tShowerView,"nonfinite_count",&m_showerView.nonfinite_count,"I");
    m_tIsoConstituent=tree("RJIsolationConstituentV1"); idBranches(m_tIsoConstituent,"candidate_id",&m_isoConstituentCandidateId); idBranches(m_tIsoConstituent,"constituent_id",&m_isoConstituentId); scalar(m_tIsoConstituent,"delta_eta",&m_isoConstituent.delta_eta,"D"); scalar(m_tIsoConstituent,"delta_phi",&m_isoConstituent.delta_phi,"D"); scalar(m_tIsoConstituent,"delta_r",&m_isoConstituent.delta_r,"D"); scalar(m_tIsoConstituent,"subsystem",&m_isoConstituent.subsystem,"I"); scalar(m_tIsoConstituent,"raw_energy",&m_isoConstituent.raw_energy,"D"); scalar(m_tIsoConstituent,"calibrated_energy",&m_isoConstituent.calibrated_energy,"D"); scalar(m_tIsoConstituent,"sub1_energy",&m_isoConstituent.sub1_energy,"D"); scalar(m_tIsoConstituent,"phosub_residual",&m_isoConstituent.phosub_residual,"D"); scalar(m_tIsoConstituent,"quality_state",&m_isoConstituent.quality_state,"I"); scalar(m_tIsoConstituent,"mask_state",&m_isoConstituent.mask_state,"I"); scalar(m_tIsoConstituent,"candidate_removal_state",&m_isoConstituent.candidate_removal_state,"I");
    scalar(m_tIsoConstituent,"native_source",&m_isoConstituent.native_source,"I"); scalar(m_tIsoConstituent,"native_key",&m_isoConstituent.native_key,"i");
    m_tIsoWitness=tree("RJIsolationWitnessV1"); idBranches(m_tIsoWitness,"candidate_id",&m_isoWitnessCandidateId); idBranches(m_tIsoWitness,"isolation_identity",&m_isoWitnessId); scalar(m_tIsoWitness,"radius",&m_isoWitness.radius,"D"); scalar(m_tIsoWitness,"subtraction_method",&m_isoWitness.subtraction_method,"I"); scalar(m_tIsoWitness,"reconstructed_or_truth",&m_isoWitness.reconstructed_or_truth,"I"); scalar(m_tIsoWitness,"cone_sum",&m_isoWitness.cone_sum,"D"); scalar(m_tIsoWitness,"threshold",&m_isoWitness.threshold,"D"); scalar(m_tIsoWitness,"sideband_threshold",&m_isoWitness.sideband_threshold,"D"); scalar(m_tIsoWitness,"pass_state",&m_isoWitness.pass_state,"I"); scalar(m_tIsoWitness,"constituent_count",&m_isoWitness.constituent_count,"I");
    m_tJet=tree("RJJetV1"); idBranches(m_tJet,"jet_id",&m_jetId); idBranches(m_tJet,"event_id",&m_jetEventId); str(m_tJet,"algorithm",&m_jet.algorithm); scalar(m_tJet,"radius",&m_jet.radius,"D"); str(m_tJet,"input_identity",&m_jet.input_identity); str(m_tJet,"subtraction_identity",&m_jet.subtraction_identity); scalar(m_tJet,"raw_pt",&m_jet.raw_pt,"D"); scalar(m_tJet,"corrected_pt",&m_jet.corrected_pt,"D"); scalar(m_tJet,"eta",&m_jet.eta,"D"); scalar(m_tJet,"phi",&m_jet.phi,"D"); scalar(m_tJet,"mass",&m_jet.mass,"D"); scalar(m_tJet,"area",&m_jet.area,"D"); unsigned64(m_tJet,"quality_bitmask",&m_jetQualityBitmask); scalar(m_tJet,"deterministic_order",&m_jet.deterministic_order,"I");
    scalar(m_tJet,"native_jet_key",&m_jet.native_jet_key,"i"); scalar(m_tJet,"native_raw_jet_key",&m_jet.native_raw_jet_key,"i");
    m_tJetConstituent=tree("RJJetConstituentV1"); idBranches(m_tJetConstituent,"jet_id",&m_jetConstituentJetId); idBranches(m_tJetConstituent,"constituent_identity",&m_jetConstituentId); scalar(m_tJetConstituent,"constituent_ordinal",&m_jetConstituent.constituent_ordinal,"I"); scalar(m_tJetConstituent,"subsystem",&m_jetConstituent.subsystem,"I"); scalar(m_tJetConstituent,"energy",&m_jetConstituent.energy,"D"); scalar(m_tJetConstituent,"eta",&m_jetConstituent.eta,"D"); scalar(m_tJetConstituent,"phi",&m_jetConstituent.phi,"D"); scalar(m_tJetConstituent,"quality_state",&m_jetConstituent.quality_state,"I");
    scalar(m_tJetConstituent,"native_source",&m_jetConstituent.native_source,"I"); scalar(m_tJetConstituent,"native_key",&m_jetConstituent.native_key,"i");
    m_tPair=tree("RJPhotonJetPairV1"); idBranches(m_tPair,"pair_id",&m_pairId); idBranches(m_tPair,"event_id",&m_pairEventId); idBranches(m_tPair,"candidate_id",&m_pairCandidateId); idBranches(m_tPair,"jet_id",&m_pairJetId); scalar(m_tPair,"delta_phi",&m_pair.delta_phi,"D"); scalar(m_tPair,"xjgamma",&m_pair.xjgamma,"D"); scalar(m_tPair,"recoil_state",&m_pair.recoil_state,"I"); scalar(m_tPair,"photon_rank",&m_pair.photon_rank,"I"); scalar(m_tPair,"jet_rank",&m_pair.jet_rank,"I"); scalar(m_tPair,"wrong_photon_class",&m_pair.wrong_photon_class,"I"); scalar(m_tPair,"wrong_recoil_class",&m_pair.wrong_recoil_class,"I");
    m_tTruthPhoton=tree("RJTruthPhotonV1"); idBranches(m_tTruthPhoton,"truth_photon_id",&m_truthPhotonId); idBranches(m_tTruthPhoton,"event_id",&m_truthPhotonEventId); scalar(m_tTruthPhoton,"pt",&m_truthPhoton.pt,"D"); scalar(m_tTruthPhoton,"eta",&m_truthPhoton.eta,"D"); scalar(m_tTruthPhoton,"phi",&m_truthPhoton.phi,"D"); scalar(m_tTruthPhoton,"prompt_class",&m_truthPhoton.prompt_class,"I"); scalar(m_tTruthPhoton,"source_role",&m_truthPhoton.source_role,"I"); scalar(m_tTruthPhoton,"g4_photon_valid",&m_truthPhoton.g4_photon_valid,"I"); scalar(m_tTruthPhoton,"hepmc_association_valid",&m_truthPhoton.hepmc_association_valid,"I"); scalar(m_tTruthPhoton,"truth_isolation_valid",&m_truthPhoton.truth_isolation_valid,"I"); scalar(m_tTruthPhoton,"analysis_signal_r03",&m_truthPhoton.analysis_signal_r03,"I"); scalar(m_tTruthPhoton,"generator_barcode",&m_truthPhoton.generator_barcode,"I"); scalar(m_tTruthPhoton,"truth_isolation_witness",&m_truthPhoton.truth_isolation_witness,"D"); scalar(m_tTruthPhoton,"truth_isolation_r03",&m_truthPhoton.truth_isolation_r03,"D"); scalar(m_tTruthPhoton,"truth_isolation_r04",&m_truthPhoton.truth_isolation_r04,"D"); scalar(m_tTruthPhoton,"reporting_guard_state",&m_truthPhoton.reporting_guard_state,"I");
    scalar(m_tTruthPhoton,"generator_occurrence_embedding_id",&m_truthPhoton.generator_occurrence_embedding_id,"I");
    scalar(m_tTruthPhoton,"sample_source_role",&m_truthPhoton.sample_source_role,"I");
    scalar(m_tTruthPhoton,"native_track_id",&m_truthPhoton.native_track_id,"I"); scalar(m_tTruthPhoton,"embedding_id",&m_truthPhoton.embedding_id,"I");
    m_tTruthPhotonMissOccurrence=tree("RJTruthPhotonMissOccurrenceV1"); idBranches(m_tTruthPhotonMissOccurrence,"occurrence_id",&m_truthPhotonMissOccurrenceId); idBranches(m_tTruthPhotonMissOccurrence,"event_id",&m_truthPhotonMissOccurrenceEventId); str(m_tTruthPhotonMissOccurrence,"trigger",&m_truthPhotonMissOccurrence.trigger); str(m_tTruthPhotonMissOccurrence,"iso_view",&m_truthPhotonMissOccurrence.iso_view); scalar(m_tTruthPhotonMissOccurrence,"centrality_bin_index",&m_truthPhotonMissOccurrence.centrality_bin_index,"I"); scalar(m_tTruthPhotonMissOccurrence,"truth_pt",&m_truthPhotonMissOccurrence.truth_pt,"D"); str(m_tTruthPhotonMissOccurrence,"decision_path",&m_truthPhotonMissOccurrence.decision_path); scalar(m_tTruthPhotonMissOccurrence,"occurrence_ordinal",&m_truthPhotonMissOccurrence.occurrence_ordinal,"I");
    m_tEmbeddedPhotonDiagnosticOccurrence=tree("RJEmbeddedPhotonDiagnosticOccurrenceV1"); idBranches(m_tEmbeddedPhotonDiagnosticOccurrence,"occurrence_id",&m_embeddedPhotonDiagnosticOccurrenceId); idBranches(m_tEmbeddedPhotonDiagnosticOccurrence,"event_id",&m_embeddedPhotonDiagnosticOccurrenceEventId); str(m_tEmbeddedPhotonDiagnosticOccurrence,"trigger",&m_embeddedPhotonDiagnosticOccurrence.trigger); str(m_tEmbeddedPhotonDiagnosticOccurrence,"occurrence_kind",&m_embeddedPhotonDiagnosticOccurrence.occurrence_kind); scalar(m_tEmbeddedPhotonDiagnosticOccurrence,"audit_bin",&m_embeddedPhotonDiagnosticOccurrence.audit_bin,"I"); scalar(m_tEmbeddedPhotonDiagnosticOccurrence,"sample_code",&m_embeddedPhotonDiagnosticOccurrence.sample_code,"I"); scalar(m_tEmbeddedPhotonDiagnosticOccurrence,"centrality_bin_index",&m_embeddedPhotonDiagnosticOccurrence.centrality_bin_index,"I"); scalar(m_tEmbeddedPhotonDiagnosticOccurrence,"photon_class",&m_embeddedPhotonDiagnosticOccurrence.photon_class,"I"); scalar(m_tEmbeddedPhotonDiagnosticOccurrence,"truth_pt",&m_embeddedPhotonDiagnosticOccurrence.truth_pt,"D"); scalar(m_tEmbeddedPhotonDiagnosticOccurrence,"truth_isolation_witness",&m_embeddedPhotonDiagnosticOccurrence.truth_isolation_witness,"D"); scalar(m_tEmbeddedPhotonDiagnosticOccurrence,"stitch_decision",&m_embeddedPhotonDiagnosticOccurrence.stitch_decision,"I"); scalar(m_tEmbeddedPhotonDiagnosticOccurrence,"stitch_photon_pt",&m_embeddedPhotonDiagnosticOccurrence.stitch_photon_pt,"D"); scalar(m_tEmbeddedPhotonDiagnosticOccurrence,"occurrence_ordinal",&m_embeddedPhotonDiagnosticOccurrence.occurrence_ordinal,"I");
    m_tPPG12DiagnosticOccurrence=tree("RJPPG12DiagnosticOccurrenceV1"); idBranches(m_tPPG12DiagnosticOccurrence,"occurrence_id",&m_ppg12DiagnosticOccurrenceId); idBranches(m_tPPG12DiagnosticOccurrence,"event_id",&m_ppg12DiagnosticOccurrenceEventId); str(m_tPPG12DiagnosticOccurrence,"trigger",&m_ppg12DiagnosticOccurrence.trigger); str(m_tPPG12DiagnosticOccurrence,"occurrence_kind",&m_ppg12DiagnosticOccurrence.occurrence_kind); scalar(m_tPPG12DiagnosticOccurrence,"stage_code",&m_ppg12DiagnosticOccurrence.stage_code,"I"); scalar(m_tPPG12DiagnosticOccurrence,"source_code",&m_ppg12DiagnosticOccurrence.source_code,"I"); scalar(m_tPPG12DiagnosticOccurrence,"sample_code",&m_ppg12DiagnosticOccurrence.sample_code,"I"); scalar(m_tPPG12DiagnosticOccurrence,"decision",&m_ppg12DiagnosticOccurrence.decision,"I"); scalar(m_tPPG12DiagnosticOccurrence,"cluster_ordinal",&m_ppg12DiagnosticOccurrence.cluster_ordinal,"I"); scalar(m_tPPG12DiagnosticOccurrence,"truth_track_id",&m_ppg12DiagnosticOccurrence.truth_track_id,"I"); scalar(m_tPPG12DiagnosticOccurrence,"generator_barcode",&m_ppg12DiagnosticOccurrence.generator_barcode,"I"); scalar(m_tPPG12DiagnosticOccurrence,"source_photon_pt",&m_ppg12DiagnosticOccurrence.source_photon_pt,"D"); scalar(m_tPPG12DiagnosticOccurrence,"flow_value",&m_ppg12DiagnosticOccurrence.flow_value,"D"); scalar(m_tPPG12DiagnosticOccurrence,"reco_photon_pt",&m_ppg12DiagnosticOccurrence.reco_photon_pt,"D"); scalar(m_tPPG12DiagnosticOccurrence,"response_photon_pt",&m_ppg12DiagnosticOccurrence.response_photon_pt,"D"); scalar(m_tPPG12DiagnosticOccurrence,"truth_photon_pt",&m_ppg12DiagnosticOccurrence.truth_photon_pt,"D"); scalar(m_tPPG12DiagnosticOccurrence,"truth_prior_weight",&m_ppg12DiagnosticOccurrence.truth_prior_weight,"D"); scalar(m_tPPG12DiagnosticOccurrence,"reco_vertex_z",&m_ppg12DiagnosticOccurrence.reco_vertex_z,"D"); scalar(m_tPPG12DiagnosticOccurrence,"hard_truth_vertex_z",&m_ppg12DiagnosticOccurrence.hard_truth_vertex_z,"D"); scalar(m_tPPG12DiagnosticOccurrence,"mb_truth_vertex_z",&m_ppg12DiagnosticOccurrence.mb_truth_vertex_z,"D"); scalar(m_tPPG12DiagnosticOccurrence,"vertex_weight",&m_ppg12DiagnosticOccurrence.vertex_weight,"D"); scalar(m_tPPG12DiagnosticOccurrence,"period_event_weight",&m_ppg12DiagnosticOccurrence.period_event_weight,"D"); scalar(m_tPPG12DiagnosticOccurrence,"cluster_energy",&m_ppg12DiagnosticOccurrence.cluster_energy,"D"); scalar(m_tPPG12DiagnosticOccurrence,"cluster_eta",&m_ppg12DiagnosticOccurrence.cluster_eta,"D"); scalar(m_tPPG12DiagnosticOccurrence,"cluster_et",&m_ppg12DiagnosticOccurrence.cluster_et,"D"); scalar(m_tPPG12DiagnosticOccurrence,"truth_energy",&m_ppg12DiagnosticOccurrence.truth_energy,"D"); scalar(m_tPPG12DiagnosticOccurrence,"truth_eta",&m_ppg12DiagnosticOccurrence.truth_eta,"D"); scalar(m_tPPG12DiagnosticOccurrence,"truth_et",&m_ppg12DiagnosticOccurrence.truth_et,"D"); scalar(m_tPPG12DiagnosticOccurrence,"energy_contribution",&m_ppg12DiagnosticOccurrence.energy_contribution,"D"); scalar(m_tPPG12DiagnosticOccurrence,"occurrence_weight",&m_ppg12DiagnosticOccurrence.occurrence_weight,"D"); unsigned64(m_tPPG12DiagnosticOccurrence,"selection_bitmask",&m_ppg12DiagnosticSelectionBitmask); scalar(m_tPPG12DiagnosticOccurrence,"occurrence_ordinal",&m_ppg12DiagnosticOccurrence.occurrence_ordinal,"I");
    m_tTruthJet=tree("RJTruthJetV1"); idBranches(m_tTruthJet,"truth_jet_id",&m_truthJetId); idBranches(m_tTruthJet,"event_id",&m_truthJetEventId); str(m_tTruthJet,"algorithm",&m_truthJet.algorithm); scalar(m_tTruthJet,"radius",&m_truthJet.radius,"D"); scalar(m_tTruthJet,"pt",&m_truthJet.pt,"D"); scalar(m_tTruthJet,"eta",&m_truthJet.eta,"D"); scalar(m_tTruthJet,"phi",&m_truthJet.phi,"D"); str(m_tTruthJet,"ownership_state",&m_truthJet.ownership_state); scalar(m_tTruthJet,"reporting_guard_state",&m_truthJet.reporting_guard_state,"I");
    scalar(m_tTruthJet,"native_jet_key",&m_truthJet.native_jet_key,"i");
    scalar(m_tTruthPhoton,"native_vertex_id",&m_truthPhoton.native_vertex_id,"I");
    m_tLink=tree("RJRecoTruthLinkV1"); idBranches(m_tLink,"link_id",&m_linkId); scalar(m_tLink,"reco_type",&m_link.reco_type,"I"); idBranches(m_tLink,"reco_id",&m_linkRecoId); scalar(m_tLink,"truth_type",&m_link.truth_type,"I"); idBranches(m_tLink,"truth_id",&m_linkTruthId); scalar(m_tLink,"match_metric",&m_link.match_metric,"D"); scalar(m_tLink,"link_class",&m_link.link_class,"I");
    m_tWeight=tree("RJWeightComponentV1"); idBranches(m_tWeight,"target_id",&m_weightTargetId); str(m_tWeight,"component_type",&m_weight.component_type); scalar(m_tWeight,"slice_weight",&m_weight.slice_weight,"D"); scalar(m_tWeight,"cross_section_weight",&m_weight.cross_section_weight,"D"); scalar(m_tWeight,"vertex_weight",&m_weight.vertex_weight,"D"); scalar(m_tWeight,"si_di_weight",&m_weight.si_di_weight,"D"); scalar(m_tWeight,"period_weight",&m_weight.period_weight,"D"); scalar(m_tWeight,"exposure_weight",&m_weight.exposure_weight,"D"); scalar(m_tWeight,"final_weight",&m_weight.final_weight,"D"); scalar(m_tWeight,"application_count",&m_weight.application_count,"I");
    scalar(m_tWeight,"target_type",&m_weight.target_type,"I");
    m_tSnapshot=tree("RJEventDisplaySnapshotV1"); idBranches(m_tSnapshot,"snapshot_id",&m_snapshotId); idBranches(m_tSnapshot,"event_id",&m_snapshotEventId); str(m_tSnapshot,"selection_reason",&m_snapshot.selection_reason); str(m_tSnapshot,"quota_class",&m_snapshot.quota_class); str(m_tSnapshot,"serialized_payload_hash",&m_snapshot.serialized_payload_hash);
  }

  bool validRecoTruthLink(const RecoTruthLinkRow& r) const
  {
    const auto rt=static_cast<RecoTruthType>(r.reco_type), tt=static_cast<RecoTruthType>(r.truth_type);
    const auto lc=static_cast<LinkClass>(r.link_class);
    const bool recoOk=(rt==RecoTruthType::PHOTON&&contains(m_candidateIds,r.reco_id))||(rt==RecoTruthType::JET&&contains(m_jetIds,r.reco_id))||(rt==RecoTruthType::NONE&&r.reco_id.isNull());
    const bool truthOk=(tt==RecoTruthType::PHOTON&&contains(m_truthPhotonIds,r.truth_id))||(tt==RecoTruthType::JET&&contains(m_truthJetIds,r.truth_id))||(tt==RecoTruthType::NONE&&r.truth_id.isNull());
    if(!recoOk||!truthOk) return false;
    if(lc==LinkClass::RECO_FAKE) return rt!=RecoTruthType::NONE&&tt==RecoTruthType::NONE;
    if(lc==LinkClass::TRUTH_MISS) return rt==RecoTruthType::NONE&&tt!=RecoTruthType::NONE;
    return rt!=RecoTruthType::NONE&&tt!=RecoTruthType::NONE;
  }

  TFile* m_file=nullptr; TDirectory* m_dir=nullptr; bool m_initialized=false,m_finished=false; WriterMode m_mode=WriterMode::SERIALIZE; Metadata m_metadata;
  std::vector<TTree*> m_trees;
  TTree *m_tSource=nullptr,*m_tEvent=nullptr,*m_tCandidate=nullptr,*m_tModel=nullptr,*m_tShower=nullptr,*m_tShowerView=nullptr,*m_tIsoConstituent=nullptr,*m_tIsoWitness=nullptr,*m_tJet=nullptr,*m_tJetConstituent=nullptr,*m_tPair=nullptr,*m_tTruthPhoton=nullptr,*m_tTruthPhotonMissOccurrence=nullptr,*m_tEmbeddedPhotonDiagnosticOccurrence=nullptr,*m_tPPG12DiagnosticOccurrence=nullptr,*m_tTruthJet=nullptr,*m_tLink=nullptr,*m_tWeight=nullptr,*m_tSnapshot=nullptr;
  SourceOccurrenceRow m_source; EventRow m_event; PhotonCandidateRow m_candidate; ModelEvaluationRow m_model; ShowerCellRow m_shower; ShowerFeatureViewRow m_showerView; IsolationConstituentRow m_isoConstituent; IsolationWitnessRow m_isoWitness; JetRow m_jet; JetConstituentRow m_jetConstituent; PhotonJetPairRow m_pair; TruthPhotonRow m_truthPhoton; TruthPhotonMissOccurrenceRow m_truthPhotonMissOccurrence; EmbeddedPhotonDiagnosticOccurrenceRow m_embeddedPhotonDiagnosticOccurrence; PPG12DiagnosticOccurrenceRow m_ppg12DiagnosticOccurrence; TruthJetRow m_truthJet; RecoTruthLinkRow m_link; WeightComponentRow m_weight; EventDisplaySnapshotRow m_snapshot;
  SerializedIdentity128 m_sourceId;
  SerializedIdentity128 m_eventId,m_eventSourceId;
  SerializedIdentity128 m_candidateId,m_candidateEventId;
  SerializedIdentity128 m_modelCandidateId,m_modelId;
  SerializedIdentity128 m_showerCandidateId;
  SerializedIdentity128 m_showerViewCandidateId,m_showerViewDefinitionId;
  SerializedIdentity128 m_isoConstituentCandidateId,m_isoConstituentId;
  SerializedIdentity128 m_isoWitnessCandidateId,m_isoWitnessId;
  SerializedIdentity128 m_jetId,m_jetEventId;
  SerializedIdentity128 m_jetConstituentJetId,m_jetConstituentId;
  SerializedIdentity128 m_pairId,m_pairEventId,m_pairCandidateId,m_pairJetId;
  SerializedIdentity128 m_truthPhotonId,m_truthPhotonEventId;
  SerializedIdentity128 m_truthPhotonMissOccurrenceId,m_truthPhotonMissOccurrenceEventId;
  SerializedIdentity128 m_embeddedPhotonDiagnosticOccurrenceId,m_embeddedPhotonDiagnosticOccurrenceEventId;
  SerializedIdentity128 m_ppg12DiagnosticOccurrenceId,m_ppg12DiagnosticOccurrenceEventId;
  SerializedIdentity128 m_truthJetId,m_truthJetEventId;
  SerializedIdentity128 m_linkId,m_linkRecoId,m_linkTruthId;
  SerializedIdentity128 m_weightTargetId;
  SerializedIdentity128 m_snapshotId,m_snapshotEventId;
  SerializedInt64 m_eventSourceFileOrdinal=-1,m_eventSourceEntryOrdinal=-1,m_eventSourceGlobalEntryOrdinal=-1;
  SerializedInt64 m_eventSequence=0,m_eventPhysicalEventSequence=0,m_eventWeightedFillCount=0,m_eventRawFillCount=0;
  std::array<SerializedInt64,5> m_eventWeightedFillCountByCode{};
  std::array<SerializedInt64,5> m_eventRawFillCountByCode{};
  SerializedUInt64 m_eventTriggerBits=0,m_eventLiveTriggerBits=0,m_eventScaledTriggerBits=0,m_eventTriggerPacketStatus=0,m_candidatePreselectionBitmask=0,m_ppg12DiagnosticSelectionBitmask=0;
  SerializedUInt64 m_eventTriggerScalerSnapshotId=0;
  SerializedUInt64 m_showerTowerKey=0,m_jetQualityBitmask=0;
  IdentityRefMap m_sourceRefs,m_eventRefs,m_candidateRefs,m_modelRefs,
      m_showerDefinitionRefs,m_jetRefs,m_truthPhotonRefs,m_truthJetRefs;
  std::uint64_t m_eventRows=0,m_candidateRows=0,m_jetRows=0,
      m_truthPhotonRows=0,m_truthJetRows=0,
      m_isoConstituentRows=0,m_isoWitnessRows=0,
      m_jetConstituentRows=0,m_pairRows=0,m_occurrenceRows=0,m_linkRows=0,
      m_snapshotRows=0;
  std::unordered_set<Identity128,IdentityHash> m_sourceIds,m_seenEventIds,
      m_eventIds,m_candidateIds,m_modelEvalIds,m_showerCellIds,m_showerViewIds,
      m_jetIds,m_pairIds,m_truthPhotonIds,m_truthPhotonMissOccurrenceIds,
      m_embeddedPhotonDiagnosticOccurrenceIds,m_ppg12DiagnosticOccurrenceIds,
      m_truthJetIds,m_linkIds,m_snapshotIds;
};
} // namespace RJReplayFoundationV1

#endif
