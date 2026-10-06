#ifndef TREEPRODUCTION_EVENTTREE_H
#define TREEPRODUCTION_EVENTTREE_H

/**
 * @file EventTree.h
 * ROOT operations on the single-eventTree product, shared by production,
 * source assembly and PhotonID. This library links ROOT only: reading,
 * checking, assembling or extending a saved tree never loads Fun4All or a
 * detector library.
 *
 * Requests and results that carry structure are JSON text (UTF-8), so
 * Python callers bind through PyROOT without an extra layer. A violated
 * contract throws std::runtime_error with the reason; outputs are written to
 * new candidate paths only and are never published by a failed call.
 */

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

class TFile;
class TObject;
class TTree;

namespace eventtree
{

/// SHA-256 of a regular file, streamed in bounded blocks (lowercase hex).
std::string fileSha256(const std::string& path);
/// SHA-256 of text, such as a canonical JSON definition.
std::string textSha256(const std::string& text);

// ============================================================================
// Schema
//
// config/eventTree_v1.yaml is the single authority. The build expands it into
// the concrete branch table compiled into this library; Schema::load refuses
// a schema file whose digest differs from that table. Rebuild after editing
// the schema.
// ============================================================================

struct Branch
{
  std::string name;
  std::string group;
  std::string type;     // storage type: int8 ... uint64, float, double, string
  std::string title;    // with its unit, as written to the TBranch
  std::string counter;  // Int_t length branch of a numeric collection
  int length = 0;       // 0 scalar, -1 collection, n > 0 fixed array
  bool reserved = false;
};

class Schema
{
 public:
  /// The compiled schema, after checking that `path` is its source file.
  static const Schema& load(const std::string& path);

  const std::vector<Branch>& branches() const { return m_branches; }
  const Branch* find(const std::string& name) const;
  const std::string& sha256() const { return m_sha256; }
  const std::string& contractSha256() const { return m_contractSha256; }

  /// This revision, or a pinned earlier revision whose complete output
  /// contract is identical. Arbitrary other digests are refused.
  bool acceptsDigest(const std::string& digest) const;

 private:
  Schema();
  std::vector<Branch> m_branches;
  std::map<std::string, std::size_t> m_index;
  std::string m_sha256;
  std::string m_contractSha256;
};

/// Photon-ID model input name -> the eventTree branch that stores it.
/// shv_ branches are read at the model's shower view; vtx_z and cent_value
/// are per event. This is the one table shared by checks and PhotonID.
const std::vector<std::pair<std::string, std::string>>& featureBranches();
std::string featureBranch(const std::string& input);

// ============================================================================
// Writer
//
// Bounded TTree writing from typed values. A collection is a C array whose
// length is an Int_t counter branch, so an empty collection costs almost
// nothing per entry. Metadata are ordinary ROOT objects, never extra rows.
// ============================================================================

class Column;

class Writer
{
 public:
  Writer(const std::string& path, const std::string& treeName, const std::string& title,
         int compression, long long entries);
  ~Writer();
  Writer(const Writer&) = delete;
  Writer& operator=(const Writer&) = delete;

  /// Rebind every column for the next batch; previous fills have finished.
  void beginBatch(long long entries);

  template <class T>
  void scalar(const std::string& name, const std::string& leaf, const T* data, const std::string& title);
  template <class T>
  void fixed(const std::string& name, const std::string& leaf, const T* data, int size, const std::string& title);
  void counter(const std::string& name, const long long* offsets);
  template <class T>
  void counted(const std::string& name, const std::string& leaf, const T* data, const long long* offsets,
               const std::string& counter, long long maxCount, const std::string& title);
  void strings(const std::string& name, std::vector<std::string> data, std::vector<long long> offsets,
               const std::string& title);

  long long fill();

  /// Text object at a path such as "meta/schema"; directories are created.
  void text(const std::string& path, const std::string& value);
  void object(const std::string& path, TObject& value);

  /// Write the single eventTree key and close; any ROOT write error throws.
  void close();

 private:
  void retitle(const std::string& name, const std::string& title);
  std::unique_ptr<TFile> m_file;
  TTree* m_tree = nullptr;
  std::vector<std::unique_ptr<Column>> m_columns;
  long long m_entries;
};

// One event's typed fields. Schema declarations determine storage width;
// integer identities never pass through floating point. Every field needs an
// explicit filling rule, including empty collections.
struct FieldSpec
{
  std::string name, type, title, counter;
  int length = 0;  // 0 scalar, -1 counted collection, positive fixed array
};

using FieldValues = std::variant<std::vector<std::int8_t>, std::vector<std::uint8_t>,
    std::vector<std::int16_t>, std::vector<std::uint16_t>, std::vector<std::int32_t>,
    std::vector<std::uint32_t>, std::vector<Long64_t>, std::vector<ULong64_t>,
    std::vector<float>, std::vector<double>, std::vector<std::string>>;

struct Field
{
  FieldSpec spec;
  FieldValues values;
  std::array<long long, 2> offsets{{0, 0}};
  bool assigned = false;
  std::string leaf;

  explicit Field(FieldSpec s);

  template <class T>
  void append(T value);
  void reset();
  void bind(Writer& writer);
};

class Entry
{
 public:
  explicit Entry(const std::vector<FieldSpec>& specs);
  void reset();
  template <class T>
  void put(const std::string& name, T value) { at(name).append(value); }
  template <class Values>
  void putAll(const std::string& name, const Values& values)
  {
    auto& field = at(name);
    field.assigned = true;
    for (const auto& value : values) field.append(value);
  }
  void empty(const std::string& name);
  std::vector<std::string> missing() const;
  /// Refuse any unfilled field, then write one entry (or book only).
  void write(Writer& writer, bool bookOnly = false);

 private:
  Field& at(const std::string& name);
  std::map<std::string, Field> m_fields;
  std::vector<std::string> m_order;
};

template <class T>
void Field::append(T value)
{
  std::visit([&](auto& data) {
    using Value = typename std::decay_t<decltype(data)>::value_type;
    if constexpr (std::is_arithmetic_v<T> && std::is_arithmetic_v<Value>)
    {
      if constexpr (std::is_integral_v<Value>)
      {
        bool invalid = false;
        if constexpr (std::is_integral_v<T>)
        {
          if constexpr (std::is_signed_v<T>)
          {
            if (value < 0)
            {
              if constexpr (!std::is_signed_v<Value>) invalid = true;
              else invalid = static_cast<std::intmax_t>(value) < std::numeric_limits<Value>::lowest();
            }
            else
              invalid = static_cast<std::uintmax_t>(value) > static_cast<std::uintmax_t>(std::numeric_limits<Value>::max());
          }
          else
            invalid = static_cast<std::uintmax_t>(value) > static_cast<std::uintmax_t>(std::numeric_limits<Value>::max());
        }
        else if constexpr (std::is_floating_point_v<T>)
        {
          const long double wide = static_cast<long double>(value);
          const long double upper = std::ldexp(1.0L, std::numeric_limits<Value>::digits);
          const long double lower = std::is_signed_v<Value> ? -upper : 0.0L;
          invalid = !std::isfinite(wide) || std::trunc(wide) != wide || wide < lower || wide >= upper;
        }
        if (invalid) throw std::runtime_error("integer overflow in " + spec.name);
      }
      data.push_back(static_cast<Value>(value));
    }
    else if constexpr (std::is_same_v<T, std::string> && std::is_same_v<Value, std::string>)
      data.push_back(value);
    else
      throw std::runtime_error("eventTree field type mismatch: " + spec.name);
  }, values);
  assigned = true;
}

// ============================================================================
// Saved products
// ============================================================================

/// Direct meta/<name> JSON objects of a saved file, as one JSON object.
std::string metadataJson(const std::string& path);

/// Structural check of one file, or of every original chunk of a consolidated
/// file: groups, strides, counters, event-local indices, rejected entries,
/// absent branches, code dictionaries, and, when `photonIdPackage` names the
/// frozen PhotonID package, the model inputs against the named shower
/// branches. Returns {"ok", "entries", "problems", "notes", ...}.
std::string checkFile(const std::string& path, const std::string& schemaPath,
                      const std::string& photonIdPackage = "");

/// The producer's own completion record against the job that requested it:
/// provenance, source interval, event accounting and saved entry order.
std::string checkProducerOutput(const std::string& path, const std::string& schemaPath,
                                const std::string& jobJson);

/// Assemble one logical source from its pieces. The request binds each piece
/// (path, sha256, piece_id, original_unit_id, event_begin, event_end), the
/// candidate output, logical_source_id, expected_events, exclusions, the
/// schema path and an optional progress file. Payloads, event order and
/// event-local references are bit-preserved; coded fields get a reversible
/// union dictionary; per-piece metadata and objects are kept. The saved
/// candidate is reread and proven before the receipt is returned.
std::string consolidate(const std::string& requestJson);

/// Reread a consolidated file against its original pieces.
std::string proveConsolidated(const std::string& path, const std::string& schemaPath);

/// Add declared collection branches to a separate candidate of a saved file.
/// Callers compute the values; this class copies every existing branch and
/// object, writes the declared metadata and proves preservation. Python-facing
/// members return void or text so that refusals reach PyROOT as exceptions.
class BranchAddition
{
 public:
  BranchAddition();
  ~BranchAddition();
  BranchAddition(const BranchAddition&) = delete;
  BranchAddition& operator=(const BranchAddition&) = delete;

  /// {"original", "candidate", "input_sha256", "schema_path",
  ///  "added": [{"name", "type": "double", "title"}], "metadata": {...}}
  void open(const std::string& requestJson);

  /// {"entries": saved entries, "next": first entry without values}
  std::string status() const;

  /// Values for the next `events` entries of one added branch: `valueCount`
  /// flattened values and `events` per-entry counts that must sum to it.
  /// Supply every added branch, then commit.
  void append(const std::string& branch, long long events, const double* values, long long valueCount,
              const long long* counts);

  /// Copy the batch's original entries with the supplied values.
  void commit(long long events);

  /// Close, write objects and metadata, reopen and prove. Returns the proof.
  std::string finish();

 private:
  struct State;
  std::unique_ptr<State> m_state;
};

/// Reopen an original and a candidate: every original branch, entry order,
/// type and carried object is unchanged; only declared branches and metadata
/// differ. {"original", "candidate", "added": [...], "metadata": {...}}
std::string provePreserved(const std::string& requestJson);

/// Resolve a published group's CURRENT pointer once, with its manifest check.
/// Returns {"directory", "manifest"}.
std::string resolveGroup(const std::string& group);

/// Publish a complete group version and advance CURRENT by compare-and-swap.
/// Every file must be validated by its caller and bound by a PASS receipt;
/// hashes are rechecked here, DATA physical events must be unique across the
/// group, and a stage-3 version must keep its base's logical-file roster.
std::string publishGroup(const std::string& requestJson);

}  // namespace eventtree

#endif  // TREEPRODUCTION_EVENTTREE_H
