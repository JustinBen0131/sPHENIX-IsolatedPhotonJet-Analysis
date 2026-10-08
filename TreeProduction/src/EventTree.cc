/**
 * @file EventTree.cc
 * ROOT operations on the single-eventTree product. See EventTree.h.
 *
 * Sections: identity and JSON, schema, writer, column access, checks, source
 * assembly, branch additions, group publication. Every reader streams entries
 * with bounded buffers; whole trees are never held in memory.
 */
#include "EventTree.h"

#include <TBranch.h>
#include <TBranchElement.h>
#include <TBufferFile.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TList.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TTree.h>

#include <nlohmann/json.hpp>

#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iterator>
#include <queue>
#include <set>
#include <sstream>

namespace eventtree
{
namespace
{
struct SchemaRow
{
  const char* name;
  const char* group;
  const char* type;
  const char* title;
  const char* unit;
  bool reserved;
};
// Build product: the schema file expanded into kSchemaRows, with its digests.
#include "EventTreeSchema.h"
}  // namespace
}  // namespace eventtree

namespace
{
using json = nlohmann::json;
namespace fs = std::filesystem;

constexpr const char* kTreeName = "eventTree";
constexpr int kCompression = 404;  // LZMA level 4, the schema's storage rule

// The earlier v1 schema file also held a conversion-only sources table. Its
// removal left this exact output contract unchanged; other digests refuse.
const std::map<std::string, std::string> kHistoricalContracts = {
    {"0ce3b08f37826282644ab6e3929582e920f53918288ad573801d2bb4cb8d4a44",
     "7129b8649b9273df2ebba079fd40a7f99a69912554bc5189881982be299f83dc"}};

// Coded fields that index a dictionary entry use -1 for "no entry".
const std::array<const char*, 6> kIndexSuffixes = {
    {"_reason", "_reason_index", "_origin_index", "_evidence_index", "_scoring_node", "_scoring_evidence"}};

void require(bool condition, const std::string& message)
{
  if (!condition) throw std::runtime_error(message);
}

bool endsWith(const std::string& text, const std::string& suffix)
{
  return text.size() >= suffix.size() && text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool indexField(const std::string& name)
{
  for (const char* suffix : kIndexSuffixes)
    if (endsWith(name, suffix)) return true;
  return false;
}

// ---------------------------------------------------------------------------
// SHA-256 (FIPS 180-4), streamed
// ---------------------------------------------------------------------------
class Sha256
{
 public:
  void update(const unsigned char* data, std::size_t size)
  {
    m_bits += static_cast<std::uint64_t>(size) * 8;
    while (size > 0)
    {
      const std::size_t take = std::min(size, m_block.size() - m_used);
      std::memcpy(m_block.data() + m_used, data, take);
      m_used += take;
      data += take;
      size -= take;
      if (m_used == m_block.size())
      {
        compress(m_block.data());
        m_used = 0;
      }
    }
  }

  std::string hex()
  {
    const std::uint64_t bits = m_bits;
    const unsigned char one = 0x80, zero = 0;
    update(&one, 1);
    while (m_used != 56) update(&zero, 1);
    unsigned char length[8];
    for (int i = 0; i < 8; ++i) length[i] = static_cast<unsigned char>(bits >> (56 - 8 * i));
    update(length, 8);
    std::ostringstream out;
    for (const std::uint32_t word : m_state) out << std::hex << std::setw(8) << std::setfill('0') << word;
    return out.str();
  }

 private:
  static std::uint32_t rotate(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

  void compress(const unsigned char* p)
  {
    static constexpr std::array<std::uint32_t, 64> k = {
        {0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
         0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
         0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
         0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
         0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
         0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
         0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
         0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2}};
    std::array<std::uint32_t, 64> w{};
    for (int i = 0; i < 16; ++i)
      w[i] = (std::uint32_t(p[4 * i]) << 24) | (std::uint32_t(p[4 * i + 1]) << 16) |
             (std::uint32_t(p[4 * i + 2]) << 8) | std::uint32_t(p[4 * i + 3]);
    for (int i = 16; i < 64; ++i)
    {
      const std::uint32_t s0 = rotate(w[i - 15], 7) ^ rotate(w[i - 15], 18) ^ (w[i - 15] >> 3);
      const std::uint32_t s1 = rotate(w[i - 2], 17) ^ rotate(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    auto [a, b, c, d, e, f, g, h] = m_state;
    for (int i = 0; i < 64; ++i)
    {
      const std::uint32_t t1 = h + (rotate(e, 6) ^ rotate(e, 11) ^ rotate(e, 25)) + ((e & f) ^ (~e & g)) + k[i] + w[i];
      const std::uint32_t t2 = (rotate(a, 2) ^ rotate(a, 13) ^ rotate(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }
    const std::array<std::uint32_t, 8> add = {{a, b, c, d, e, f, g, h}};
    for (int i = 0; i < 8; ++i) m_state[i] += add[i];
  }

  std::array<std::uint32_t, 8> m_state = {
      {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19}};
  std::array<unsigned char, 64> m_block{};
  std::size_t m_used = 0;
  std::uint64_t m_bits = 0;
};


// ---------------------------------------------------------------------------
// JSON
// ---------------------------------------------------------------------------
json parse(const std::string& text, const std::string& what)
{
  try
  {
    return json::parse(text);
  }
  catch (const json::exception& error)
  {
    throw std::runtime_error("invalid JSON in " + what + ": " + error.what());
  }
}

void requireFinite(const json& value, const std::string& where)
{
  if (value.is_number_float())
    require(std::isfinite(value.get<double>()), "non-finite JSON number in " + where);
  else if (value.is_structured())
    for (const auto& item : value) requireFinite(item, where);
}

// Sorted keys, compact separators, ASCII escapes: identical to the JSON the
// Python tools write with sort_keys, (',', ':') separators and allow_nan=False.
std::string dumps(const json& value)
{
  requireFinite(value, "canonical JSON");
  return value.dump(-1, ' ', true);
}

std::string utcNow()
{
  const auto now = std::time(nullptr);
  std::ostringstream out;
  out << std::put_time(std::gmtime(&now), "%Y-%m-%dT%H:%M:%SZ");
  return out.str();
}

json readJsonFile(const fs::path& path)
{
  std::ifstream in(path, std::ios::binary);
  require(static_cast<bool>(in), "cannot read " + path.string());
  std::ostringstream text;
  text << in.rdbuf();
  return parse(text.str(), path.string());
}

void writeFileAtomically(const fs::path& path, const std::string& text)
{
  const fs::path temporary = path.string() + ".tmp";
  {
    std::ofstream out(temporary, std::ios::trunc | std::ios::binary);
    out << text;
    out.flush();
    require(static_cast<bool>(out), "cannot write " + temporary.string());
  }
  fs::rename(temporary, path);
}

// Progress for long operations, replaced atomically (never appended).
class Progress
{
 public:
  Progress(std::string path, std::string phase, std::string identity, long long total)
      : m_path(std::move(path)), m_phase(std::move(phase)), m_identity(std::move(identity)), m_total(total)
  {
  }

  void update(long long processed, bool force = false)
  {
    if (m_path.empty()) return;
    const auto now = std::chrono::steady_clock::now();
    if (!force && processed % 4096 != 0 && std::chrono::duration<double>(now - m_written).count() < 5.) return;
    const json record = {{"schema", "EventTreeProgressV1"},
                         {"phase", m_phase},
                         {"logical_source_id", m_identity},
                         {"processed", processed},
                         {"total", m_total},
                         {"remaining", m_total - processed},
                         {"elapsed_seconds", std::chrono::duration<double>(now - m_started).count()},
                         {"updated_at", utcNow()}};
    writeFileAtomically(m_path, record.dump() + "\n");
    m_written = now;
  }

 private:
  std::string m_path, m_phase, m_identity;
  long long m_total;
  std::chrono::steady_clock::time_point m_started = std::chrono::steady_clock::now(), m_written = m_started;
};

// ---------------------------------------------------------------------------
// ROOT file helpers
// ---------------------------------------------------------------------------
class NoHistogramDirectory
{
 public:
  NoHistogramDirectory() : m_status(TH1::AddDirectoryStatus()) { TH1::AddDirectory(false); }
  ~NoHistogramDirectory() { TH1::AddDirectory(m_status); }

 private:
  bool m_status;
};

std::unique_ptr<TFile> openRead(const std::string& path)
{
  require(fs::is_regular_file(path), "missing eventTree file: " + path);
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  require(file && !file->IsZombie(), "cannot open " + path);
  return file;
}

std::unique_ptr<TFile> createFile(const std::string& path)
{
  require(!fs::exists(path), "output exists: " + path);
  if (fs::path(path).has_parent_path()) fs::create_directories(fs::path(path).parent_path());
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "CREATE"));
  require(file && !file->IsZombie(), "cannot create " + path);
  file->SetCompressionSettings(kCompression);
  return file;
}

TTree* eventTree(TFile& file)
{
  auto* tree = dynamic_cast<TTree*>(file.Get(kTreeName));
  require(tree != nullptr, std::string(file.GetName()) + " has no eventTree");
  return tree;
}

// The highest cycle of a name. TDirectory::Get(name) returns the first listed
// cycle, and a file updated by another writer can list an older cycle first.
TKey* latestKey(TDirectory* directory, const std::string& name)
{
  TKey* best = nullptr;
  for (TObject* item : *directory->GetListOfKeys())
  {
    auto* key = static_cast<TKey*>(item);
    if (name == key->GetName() && (!best || key->GetCycle() > best->GetCycle())) best = key;
  }
  return best;
}

// The highest cycle of an object path such as "hist/SIM/h_centrality".
TKey* latestKeyAt(TFile& file, const std::string& path)
{
  const auto separator = path.find_last_of('/');
  TDirectory* directory = separator == std::string::npos ? &file : file.GetDirectory(path.substr(0, separator).c_str());
  return directory ? latestKey(directory, separator == std::string::npos ? path : path.substr(separator + 1)) : nullptr;
}

// Highest cycle of every key below `directory`, with its full path.
void visitKeys(TDirectory* directory, const std::string& prefix,
               const std::function<void(const std::string&, TKey*)>& visit)
{
  std::set<std::string> seen;
  for (TObject* item : *directory->GetListOfKeys())
  {
    auto* key = static_cast<TKey*>(item);
    if (!seen.insert(key->GetName()).second) continue;
    key = latestKey(directory, key->GetName());
    const std::string path = prefix + key->GetName();
    TClass* type = TClass::GetClass(key->GetClassName());
    if (type && type->InheritsFrom(TDirectory::Class()))
    {
      auto* sub = directory->GetDirectory(key->GetName());
      require(sub != nullptr, "cannot enter directory " + path);
      visitKeys(sub, path + "/", visit);
    }
    else
      visit(path, key);
  }
}

bool isTree(TKey* key)
{
  TClass* type = TClass::GetClass(key->GetClassName());
  return type && type->InheritsFrom(TTree::Class());
}

// Every TTree key and cycle: a valid product holds exactly one, eventTree.
std::vector<std::string> treeKeys(TFile& file)
{
  std::vector<std::string> keys;
  std::function<void(TDirectory*, const std::string&)> walk = [&](TDirectory* directory, const std::string& prefix)
  {
    for (TObject* item : *directory->GetListOfKeys())
    {
      auto* key = static_cast<TKey*>(item);
      TClass* type = TClass::GetClass(key->GetClassName());
      if (type && type->InheritsFrom(TDirectory::Class()))
      {
        if (auto* sub = directory->GetDirectory(key->GetName())) walk(sub, prefix + key->GetName() + "/");
      }
      else if (type && type->InheritsFrom(TTree::Class()))
        keys.push_back(prefix + key->GetName() + ";" + std::to_string(key->GetCycle()));
    }
  };
  walk(&file, "");
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  return keys;
}

TDirectory* makeDirectory(TFile& file, const std::string& directory)
{
  if (directory.empty()) return &file;
  TDirectory* target = file.GetDirectory(directory.c_str());
  if (!target)
  {
    file.mkdir(directory.c_str(), "", true);
    target = file.GetDirectory(directory.c_str());
  }
  require(target != nullptr, "cannot create directory " + directory);
  return target;
}

void writeObject(TFile& file, const std::string& path, TObject& value)
{
  require(!path.empty() && path != kTreeName && path.find("..") == std::string::npos, "invalid object path: " + path);
  const auto separator = path.find_last_of('/');
  TDirectory* target = makeDirectory(file, separator == std::string::npos ? "" : path.substr(0, separator));
  const std::string name = separator == std::string::npos ? path : path.substr(separator + 1);
  target->cd();
  require(value.Write(name.c_str(), TObject::kOverwrite) > 0, "cannot write object " + path);
  file.cd();
}

void writeText(TFile& file, const std::string& path, const std::string& value)
{
  TObjString text(value.c_str());
  writeObject(file, path, text);
}

// The streamed bytes of an object: the uncompressed payload of its key.
std::string streamed(TObject& object)
{
  TBufferFile buffer(TBuffer::kWrite);
  object.Streamer(buffer);
  return std::string(buffer.Buffer(), buffer.Length());
}

std::unique_ptr<TObject> readObject(TKey* key)
{
  NoHistogramDirectory guard;
  std::unique_ptr<TObject> object(key->ReadObj());
  require(object != nullptr, std::string("cannot read object ") + key->GetName());
  return object;
}

json readMetadata(TFile& file)
{
  json result = json::object();
  TDirectory* meta = file.GetDirectory("meta");
  if (!meta) return result;
  std::set<std::string> seen;
  for (TObject* item : *meta->GetListOfKeys())
  {
    const std::string name = static_cast<TKey*>(item)->GetName();
    if (!seen.insert(name).second) continue;
    TKey* key = latestKey(meta, name);
    if (std::string(key->GetClassName()) != "TObjString") continue;
    std::unique_ptr<TObjString> text(dynamic_cast<TObjString*>(key->ReadObj()));
    require(text != nullptr, "cannot read meta/" + name);
    result[name] = parse(text->GetString().Data(), "meta/" + name);
  }
  return result;
}

std::map<std::string, std::string> keyValues(TFile& file, const std::string& name)
{
  TKey* key = latestKey(&file, name);
  std::unique_ptr<TObjString> text(key ? dynamic_cast<TObjString*>(key->ReadObj()) : nullptr);
  require(text != nullptr, std::string(file.GetName()) + " has no " + name + " record");
  std::map<std::string, std::string> values;
  std::istringstream lines(text->GetString().Data());
  for (std::string line; std::getline(lines, line);)
  {
    const auto equal = line.find('=');
    if (equal == std::string::npos) continue;
    require(values.emplace(line.substr(0, equal), line.substr(equal + 1)).second,
            "duplicate " + name + " key: " + line.substr(0, equal));
  }
  return values;
}

// A regular file's identity, to detect replacement or modification.
std::array<long long, 4> fileIdentity(const std::string& path)
{
  struct stat info{};
  require(::stat(path.c_str(), &info) == 0, "cannot stat " + path);
#ifdef __APPLE__
  const long long nanoseconds =
      static_cast<long long>(info.st_mtimespec.tv_sec) * 1000000000LL + info.st_mtimespec.tv_nsec;
#else
  const long long nanoseconds = static_cast<long long>(info.st_mtim.tv_sec) * 1000000000LL + info.st_mtim.tv_nsec;
#endif
  return {{static_cast<long long>(info.st_dev), static_cast<long long>(info.st_ino),
           static_cast<long long>(info.st_size), nanoseconds}};
}

void fsyncPath(const std::string& path)
{
  const int fd = ::open(path.c_str(), O_RDONLY);
  require(fd >= 0, "cannot open for fsync: " + path);
  const int status = ::fsync(fd);
  ::close(fd);
  require(status == 0, "fsync failed: " + path);
}
}  // namespace

// ===========================================================================
// Identity and schema
// ===========================================================================
namespace eventtree
{
std::string fileSha256(const std::string& path)
{
  require(fs::is_regular_file(path), "not a regular file: " + path);
  std::ifstream in(path, std::ios::binary);
  require(static_cast<bool>(in), "cannot read " + path);
  Sha256 digest;
  std::vector<char> block(8 << 20);
  while (in)
  {
    in.read(block.data(), static_cast<std::streamsize>(block.size()));
    digest.update(reinterpret_cast<const unsigned char*>(block.data()), static_cast<std::size_t>(in.gcount()));
  }
  require(in.eof(), "read error: " + path);
  return digest.hex();
}

Schema::Schema() : m_sha256(EVENTTREE_SCHEMA_SHA256), m_contractSha256(EVENTTREE_CONTRACT_SHA256)
{
  static const std::array<const char*, 3> radii = {{"02", "03", "04"}};
  std::set<std::string> jetCounts;
  for (const char* prefix : {"jet", "jetunsub", "truthjet"})
    for (const char* radius : radii) jetCounts.insert(std::string(prefix) + radius);
  for (const SchemaRow& row : kSchemaRows)
  {
    Branch branch;
    branch.name = row.name;
    branch.group = row.group;
    branch.title = row.title;
    branch.reserved = row.reserved;
    const std::string unit = row.unit, type = row.type;
    if (!unit.empty() && branch.title.find('[') == std::string::npos) branch.title += " [" + unit + "]";
    if (type.rfind("vector<", 0) == 0 && endsWith(type, ">"))
    {
      branch.type = type.substr(7, type.size() - 8);
      branch.length = -1;
    }
    else if (const auto bracket = type.find('['); bracket != std::string::npos)
    {
      branch.type = type.substr(0, bracket);
      branch.length = std::stoi(type.substr(bracket + 1));
      require(branch.length > 0, "invalid fixed length: " + branch.name);
    }
    else
      branch.type = type;
    // Every numeric collection has its own counter; the jet collections'
    // counters are named for the collection, anchored to pT.
    if (branch.length == -1 && branch.type != "string")
    {
      branch.counter = "n_" + branch.name;
      if (endsWith(branch.name, "_pt") && jetCounts.count(branch.name.substr(0, branch.name.size() - 3)))
        branch.counter = "n_" + branch.name.substr(0, branch.name.size() - 3);
    }
    require(m_index.emplace(branch.name, m_branches.size()).second, "duplicate schema branch: " + branch.name);
    m_branches.push_back(std::move(branch));
  }
}

std::string textSha256(const std::string& text)
{
  Sha256 digest;
  digest.update(reinterpret_cast<const unsigned char*>(text.data()), text.size());
  return digest.hex();
}

const Schema& Schema::load(const std::string& path)
{
  static const Schema compiled;
  const std::string digest = fileSha256(path);
  require(digest == compiled.m_sha256, "schema " + path + " (sha256 " + digest +
                                           ") is not the schema built into this library (" + compiled.m_sha256 +
                                           "); rebuild and reinstall TreeProduction");
  return compiled;
}

const Branch* Schema::find(const std::string& name) const
{
  const auto it = m_index.find(name);
  return it == m_index.end() ? nullptr : &m_branches[it->second];
}

bool Schema::acceptsDigest(const std::string& digest) const
{
  if (digest == m_sha256) return true;
  const auto it = kHistoricalContracts.find(digest);
  return it != kHistoricalContracts.end() && it->second == m_contractSha256;
}

const std::vector<std::pair<std::string, std::string>>& featureBranches()
{
  static const std::vector<std::pair<std::string, std::string>> table = {{"cluster_Et", "pho_et"},
                                                                         {"cluster_Eta", "pho_eta"},
                                                                         {"vertexz", "vtx_z"},
                                                                         {"centrality", "cent_value"},
                                                                         {"cluster_weta_cogx", "shv_weta_cogx"},
                                                                         {"cluster_wphi_cogx", "shv_wphi_cogx"},
                                                                         {"cluster_weta33_cogx", "shv_weta33_cogx"},
                                                                         {"cluster_wphi33_cogx", "shv_wphi33_cogx"},
                                                                         {"e11_over_e33", "shv_e11_over_e33"},
                                                                         {"e32_over_e35", "shv_e32_over_e35"},
                                                                         {"cluster_et1", "shv_et1"},
                                                                         {"cluster_et2", "shv_et2"},
                                                                         {"cluster_et3", "shv_et3"},
                                                                         {"cluster_et4", "shv_et4"}};
  return table;
}

std::string featureBranch(const std::string& input)
{
  for (const auto& item : featureBranches())
    if (item.first == input) return item.second;
  throw std::runtime_error("no eventTree branch stores model input " + input);
}

// ===========================================================================
// Writer
// ===========================================================================
class Column
{
 public:
  virtual ~Column() = default;
  virtual void load(long long entry) = 0;
};

namespace
{
template <typename T>
struct ScalarColumn : Column
{
  const T* data;
  T value{};
  ScalarColumn(TTree* tree, const std::string& name, const std::string& leaf, const T* d) : data(d)
  {
    if (auto* branch = tree->GetBranch(name.c_str()))
      branch->SetAddress(&value);
    else
      tree->Branch(name.c_str(), &value, (name + "/" + leaf).c_str());
  }
  void load(long long entry) override { value = data[entry]; }
};

template <typename T>
struct FixedColumn : Column
{
  const T* data;
  int size;
  std::vector<T> value;
  FixedColumn(TTree* tree, const std::string& name, const std::string& leaf, const T* d, int n)
      : data(d), size(n), value(n)
  {
    if (auto* branch = tree->GetBranch(name.c_str()))
      branch->SetAddress(value.data());
    else
      tree->Branch(name.c_str(), value.data(), (name + "[" + std::to_string(n) + "]/" + leaf).c_str());
  }
  void load(long long entry) override { std::memcpy(value.data(), data + entry * size, size * sizeof(T)); }
};

struct CounterColumn : Column
{
  const long long* offsets;
  Int_t value = 0;
  CounterColumn(TTree* tree, const std::string& name, const long long* o) : offsets(o)
  {
    if (auto* branch = tree->GetBranch(name.c_str()))
      branch->SetAddress(&value);
    else
      tree->Branch(name.c_str(), &value, (name + "/I").c_str());
  }
  void load(long long entry) override { value = static_cast<Int_t>(offsets[entry + 1] - offsets[entry]); }
};

template <typename T>
struct CountedColumn : Column
{
  const T* data;
  const long long* offsets;
  std::vector<T> buffer;
  CountedColumn(TTree* tree, const std::string& name, const std::string& leaf, const T* d, const long long* o,
                const std::string& counter, long long maxCount)
      : data(d), offsets(o), buffer(maxCount > 0 ? maxCount : 1)
  {
    if (auto* branch = tree->GetBranch(name.c_str()))
      branch->SetAddress(buffer.data());
    else
      tree->Branch(name.c_str(), buffer.data(), (name + "[" + counter + "]/" + leaf).c_str());
  }
  void load(long long entry) override { std::copy(data + offsets[entry], data + offsets[entry + 1], buffer.begin()); }
};

struct StringColumn : Column
{
  std::vector<std::string> data;
  std::vector<long long> offsets;
  std::vector<std::string> value;
  std::vector<std::string>* address = &value;
  StringColumn(TTree* tree, const std::string& name, std::vector<std::string> d, std::vector<long long> o)
      : data(std::move(d)), offsets(std::move(o))
  {
    if (auto* branch = tree->GetBranch(name.c_str()))
      branch->SetAddress(&address);
    else
      tree->Branch(name.c_str(), &address);
  }
  void load(long long entry) override
  {
    value.assign(data.begin() + offsets[entry], data.begin() + offsets[entry + 1]);
  }
};
}  // namespace

Writer::Writer(const std::string& path, const std::string& treeName, const std::string& title, int compression,
               long long entries)
    : m_entries(entries)
{
  m_file.reset(TFile::Open(path.c_str(), "CREATE"));
  require(m_file && !m_file->IsZombie(), "cannot create " + path);
  m_file->SetCompressionSettings(compression);
  m_tree = new TTree(treeName.c_str(), title.c_str());
  m_tree->SetAutoFlush(-16 * 1024 * 1024);
}

Writer::~Writer() = default;

void Writer::beginBatch(long long entries)
{
  m_tree->ResetBranchAddresses();
  m_columns.clear();
  m_entries = entries;
}

template <class T>
void Writer::scalar(const std::string& name, const std::string& leaf, const T* data, const std::string& title)
{
  m_columns.emplace_back(new ScalarColumn<T>(m_tree, name, leaf, data));
  retitle(name, title);
}

template <class T>
void Writer::fixed(const std::string& name, const std::string& leaf, const T* data, int size, const std::string& title)
{
  m_columns.emplace_back(new FixedColumn<T>(m_tree, name, leaf, data, size));
  retitle(name, title);
}

void Writer::counter(const std::string& name, const long long* offsets)
{
  m_columns.emplace_back(new CounterColumn(m_tree, name, offsets));
}

template <class T>
void Writer::counted(const std::string& name, const std::string& leaf, const T* data, const long long* offsets,
                     const std::string& counter, long long maxCount, const std::string& title)
{
  m_columns.emplace_back(new CountedColumn<T>(m_tree, name, leaf, data, offsets, counter, maxCount));
  retitle(name, title);
}

#define EVENTTREE_WRITER_TYPE(T)                                                                             \
  template void Writer::scalar<T>(const std::string&, const std::string&, const T*, const std::string&);     \
  template void Writer::fixed<T>(const std::string&, const std::string&, const T*, int, const std::string&); \
  template void Writer::counted<T>(const std::string&, const std::string&, const T*, const long long*,       \
                                   const std::string&, long long, const std::string&);
EVENTTREE_WRITER_TYPE(std::int8_t)
EVENTTREE_WRITER_TYPE(std::uint8_t)
EVENTTREE_WRITER_TYPE(std::int16_t)
EVENTTREE_WRITER_TYPE(std::uint16_t)
EVENTTREE_WRITER_TYPE(std::int32_t)
EVENTTREE_WRITER_TYPE(std::uint32_t)
EVENTTREE_WRITER_TYPE(Long64_t)
EVENTTREE_WRITER_TYPE(ULong64_t)
EVENTTREE_WRITER_TYPE(float)
EVENTTREE_WRITER_TYPE(double)
#undef EVENTTREE_WRITER_TYPE

void Writer::strings(const std::string& name, std::vector<std::string> data, std::vector<long long> offsets,
                     const std::string& title)
{
  m_columns.emplace_back(new StringColumn(m_tree, name, std::move(data), std::move(offsets)));
  retitle(name, title);
}

long long Writer::fill()
{
  for (long long entry = 0; entry < m_entries; ++entry)
  {
    for (auto& column : m_columns) column->load(entry);
    require(m_tree->Fill() >= 0, "TTree::Fill failed");
  }
  return m_entries;
}

void Writer::text(const std::string& path, const std::string& value) { writeText(*m_file, path, value); }

void Writer::object(const std::string& path, TObject& value) { writeObject(*m_file, path, value); }

void Writer::close()
{
  m_file->cd();
  // kOverwrite replaces an autosave cycle: the file holds exactly one eventTree key.
  require(m_tree->Write("", TObject::kOverwrite) > 0, "TTree::Write failed");
  require(!m_file->TestBit(TFile::kWriteError), "ROOT file write error");
  m_file->Close();
  require(!m_file->TestBit(TFile::kWriteError), "ROOT file close error");
}

void Writer::retitle(const std::string& name, const std::string& title)
{
  if (!title.empty()) m_tree->GetBranch(name.c_str())->SetTitle(title.c_str());
}

Field::Field(FieldSpec s) : spec(std::move(s))
{
  const std::string& type = spec.type;
  if (type == "int8")
  {
    values = std::vector<std::int8_t>();
    leaf = "B";
  }
  else if (type == "uint8")
  {
    values = std::vector<std::uint8_t>();
    leaf = "b";
  }
  else if (type == "int16")
  {
    values = std::vector<std::int16_t>();
    leaf = "S";
  }
  else if (type == "uint16")
  {
    values = std::vector<std::uint16_t>();
    leaf = "s";
  }
  else if (type == "int32")
  {
    values = std::vector<std::int32_t>();
    leaf = "I";
  }
  else if (type == "uint32")
  {
    values = std::vector<std::uint32_t>();
    leaf = "i";
  }
  else if (type == "int64")
  {
    values = std::vector<Long64_t>();
    leaf = "L";
  }
  else if (type == "uint64")
  {
    values = std::vector<ULong64_t>();
    leaf = "l";
  }
  else if (type == "float")
  {
    values = std::vector<float>();
    leaf = "F";
  }
  else if (type == "double")
  {
    values = std::vector<double>();
    leaf = "D";
  }
  else if (type == "string" && spec.length == -1)
    values = std::vector<std::string>();
  else
    throw std::runtime_error("unsupported eventTree field type: " + type);
}

void Field::reset()
{
  std::visit([](auto& data) { data.clear(); }, values);
  assigned = false;
}

void Field::bind(Writer& writer)
{
  std::visit(
      [&](auto& data)
      {
        using Value = typename std::decay_t<decltype(data)>::value_type;
        if (spec.length >= 0 && data.size() != static_cast<std::size_t>(std::max(1, spec.length)))
          throw std::runtime_error("scalar/fixed field length differs: " + spec.name);
        offsets[1] = static_cast<long long>(data.size());
        if constexpr (std::is_same_v<Value, std::string>)
          writer.strings(spec.name, data, {0, offsets[1]}, spec.title);
        else
        {
          static const Value empty{};
          const Value* pointer = data.empty() ? &empty : data.data();
          if (spec.length == -1)
          {
            require(data.size() <= static_cast<std::size_t>(std::numeric_limits<Int_t>::max()),
                    "collection exceeds ROOT counter range: " + spec.name);
            writer.counter(spec.counter, offsets.data());
            writer.counted(spec.name, leaf, pointer, offsets.data(), spec.counter, offsets[1], spec.title);
          }
          else if (spec.length == 0)
            writer.scalar(spec.name, leaf, pointer, spec.title);
          else
            writer.fixed(spec.name, leaf, pointer, spec.length, spec.title);
        }
      },
      values);
}

Entry::Entry(const std::vector<FieldSpec>& specs)
{
  for (const auto& spec : specs)
  {
    require(m_fields.emplace(spec.name, Field(spec)).second, "duplicate eventTree field: " + spec.name);
    m_order.push_back(spec.name);
  }
}

void Entry::reset()
{
  for (auto& item : m_fields) item.second.reset();
}

Field& Entry::at(const std::string& name)
{
  const auto it = m_fields.find(name);
  require(it != m_fields.end(), "not an eventTree field: " + name);
  return it->second;
}

void Entry::empty(const std::string& name)
{
  auto& field = at(name);
  require(field.spec.length == -1, "only vector fields may be declared empty: " + name);
  field.assigned = true;
}

std::vector<std::string> Entry::missing() const
{
  std::vector<std::string> names;
  for (const auto& item : m_fields)
    if (!item.second.assigned) names.push_back(item.first);
  return names;
}

// Default buffers never masquerade as a capture: every field needs a rule.
void Entry::write(Writer& writer, bool bookOnly)
{
  const auto gaps = missing();
  require(gaps.empty(), gaps.empty() ? "" : "unfilled eventTree field: " + gaps.front());
  writer.beginBatch(bookOnly ? 0 : 1);
  for (const auto& name : m_order) m_fields.at(name).bind(writer);
  writer.fill();
}
}  // namespace eventtree

// ===========================================================================
// Column access: every branch of a saved tree bound to a typed buffer
// ===========================================================================
namespace
{
struct Buffer
{
  enum class Kind
  {
    Scalar,
    Fixed,
    Counted,
    Strings
  };
  std::string name;
  Kind kind = Kind::Scalar;
  std::string type;  // ROOT leaf type name
  std::size_t width = 0;
  int length = 1;  // fixed length
  std::string counterName;
  const Buffer* counter = nullptr;
  std::vector<std::uint64_t> storage;  // 8-byte aligned element bytes
  std::size_t capacity = 0;
  std::unique_ptr<std::vector<std::string>> strings;
  std::vector<std::string>* address = nullptr;
  TBranch* branch = nullptr;  // in the tree being read

  unsigned char* bytes() { return reinterpret_cast<unsigned char*>(storage.data()); }
  const unsigned char* bytes() const { return reinterpret_cast<const unsigned char*>(storage.data()); }
  std::size_t count() const
  {
    if (kind == Kind::Scalar) return 1;
    if (kind == Kind::Fixed) return static_cast<std::size_t>(length);
    if (kind == Kind::Strings) return strings->size();
    const Int_t value = *reinterpret_cast<const Int_t*>(counter->bytes());
    require(value >= 0, "negative counter " + counterName);
    return static_cast<std::size_t>(value);
  }
  void reserve(std::size_t elements)
  {
    capacity = std::max<std::size_t>(elements, 1);
    storage.assign((capacity * width + 7) / 8, 0);
  }
  std::int32_t* int32s()
  {
    require(type == "Int_t", name + " is not an int32 branch");
    return reinterpret_cast<std::int32_t*>(bytes());
  }
};

std::size_t typeWidth(const std::string& type)
{
  static const std::map<std::string, std::size_t> widths = {
      {"Char_t", 1}, {"UChar_t", 1}, {"Bool_t", 1},   {"Short_t", 2},   {"UShort_t", 2}, {"Int_t", 4},
      {"UInt_t", 4}, {"Float_t", 4}, {"Long64_t", 8}, {"ULong64_t", 8}, {"Double_t", 8}};
  const auto it = widths.find(type);
  require(it != widths.end(), "unsupported leaf type " + type);
  return it->second;
}

class Columns
{
 public:
  explicit Columns(TTree* tree) : m_tree(tree)
  {
    m_tree->SetCacheSize(32 << 20);
    m_tree->AddBranchToCache("*", true);
    for (TObject* item : *tree->GetListOfBranches())
    {
      auto* branch = static_cast<TBranch*>(item);
      auto buffer = std::make_unique<Buffer>();
      buffer->name = branch->GetName();
      buffer->branch = branch;
      if (auto* element = dynamic_cast<TBranchElement*>(branch))
      {
        require(std::string(element->GetClassName()) == "vector<string>",
                "unsupported object branch " + buffer->name + " (" + element->GetClassName() + ")");
        buffer->kind = Buffer::Kind::Strings;
        buffer->strings = std::make_unique<std::vector<std::string>>();
        buffer->address = buffer->strings.get();
      }
      else
      {
        require(branch->GetListOfLeaves()->GetEntries() == 1, "multi-leaf branch " + buffer->name);
        auto* leaf = static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0));
        buffer->type = leaf->GetTypeName();
        buffer->width = typeWidth(buffer->type);
        if (TLeaf* count = leaf->GetLeafCount())
        {
          buffer->kind = Buffer::Kind::Counted;
          buffer->counterName = count->GetName();
          buffer->reserve(static_cast<std::size_t>(std::max(1, count->GetMaximum())));
        }
        else
        {
          buffer->length = leaf->GetLenStatic();
          buffer->kind = buffer->length > 1 ? Buffer::Kind::Fixed : Buffer::Kind::Scalar;
          buffer->reserve(static_cast<std::size_t>(buffer->length));
        }
      }
      m_index[buffer->name] = m_buffers.size();
      m_buffers.push_back(std::move(buffer));
    }
    for (auto& buffer : m_buffers)
      if (buffer->kind == Buffer::Kind::Counted)
      {
        Buffer* counter = find(buffer->counterName);
        require(counter && counter->kind == Buffer::Kind::Scalar && counter->type == "Int_t",
                buffer->name + " has no Int_t counter " + buffer->counterName);
        buffer->counter = counter;
      }
    bind(m_tree);
  }

  Buffer* find(const std::string& name)
  {
    const auto it = m_index.find(name);
    return it == m_index.end() ? nullptr : m_buffers[it->second].get();
  }
  Buffer& at(const std::string& name)
  {
    Buffer* buffer = find(name);
    require(buffer != nullptr, std::string(m_tree->GetCurrentFile() ? m_tree->GetCurrentFile()->GetName() : "tree") +
                                   " has no branch " + name);
    return *buffer;
  }
  const std::vector<std::unique_ptr<Buffer>>& all() const { return m_buffers; }
  std::map<std::string, std::string> layout() const;
  long long entries() const { return m_tree->GetEntries(); }
  TTree* tree() const { return m_tree; }

  // Also point an output tree's matching branches at these buffers.
  void attach(TTree* output)
  {
    m_outputs.push_back(output);
    bind(output);
  }

  void read(long long entry)
  {
    // Counters precede their collections; grow a buffer before reading into it.
    for (auto& buffer : m_buffers)
      if (buffer->kind == Buffer::Kind::Scalar || buffer->kind == Buffer::Kind::Fixed) get(*buffer, entry);
    bool grown = false;
    for (auto& buffer : m_buffers)
      if (buffer->kind == Buffer::Kind::Counted && buffer->count() > buffer->capacity)
      {
        buffer->reserve(buffer->count());
        grown = true;
      }
    if (grown)
    {
      bind(m_tree);
      for (TTree* output : m_outputs) bind(output);
    }
    for (auto& buffer : m_buffers)
      if (buffer->kind == Buffer::Kind::Counted || buffer->kind == Buffer::Kind::Strings) get(*buffer, entry);
  }

 private:
  void get(Buffer& buffer, long long entry)
  {
    require(buffer.branch->GetEntry(entry) >= 0, "cannot read " + buffer.name + " entry " + std::to_string(entry));
  }

  void bind(TTree* tree)
  {
    for (auto& buffer : m_buffers)
    {
      TBranch* branch = tree->GetBranch(buffer->name.c_str());
      if (!branch) continue;  // outputs may omit nothing but may add branches
      if (buffer->kind == Buffer::Kind::Strings)
        tree->SetBranchAddress(buffer->name.c_str(), &buffer->address);
      else
        branch->SetAddress(buffer->bytes());
    }
  }

  TTree* m_tree;
  std::vector<std::unique_ptr<Buffer>> m_buffers;
  std::map<std::string, std::size_t> m_index;
  std::vector<TTree*> m_outputs;
};

// Exact equality of one buffer's current entry: counts, then payload bytes
// (NaN payloads and signed zeros included).
bool sameValues(const Buffer& a, const Buffer& b)
{
  if (a.kind != b.kind || a.type != b.type || a.length != b.length) return false;
  if (a.kind == Buffer::Kind::Strings) return *a.strings == *b.strings;
  const std::size_t n = a.count();
  return n == b.count() && std::memcmp(a.bytes(), b.bytes(), n * a.width) == 0;
}

std::string branchSignature(const Buffer& buffer)
{
  switch (buffer.kind)
  {
    case Buffer::Kind::Strings:
      return "vector<string>";
    case Buffer::Kind::Counted:
      return buffer.type + "[]";
    case Buffer::Kind::Fixed:
      return buffer.type + "[" + std::to_string(buffer.length) + "]";
    default:
      return buffer.type;
  }
}

std::map<std::string, std::string> Columns::layout() const
{
  std::map<std::string, std::string> result;
  for (const auto& buffer : m_buffers)
    result[buffer->name] =
        branchSignature(*buffer) + (buffer->counterName.empty() ? "" : " counted by " + buffer->counterName);
  return result;
}

// Map coded values between a chunk's dictionaries and the union (or back).
void remap(Columns& columns, const std::map<std::string, std::vector<int>>& maps)
{
  for (const auto& [name, mapping] : maps)
  {
    Buffer& buffer = columns.at(name);
    std::int32_t* values = buffer.int32s();
    const long long lowest = indexField(name) ? -1 : 0;
    for (std::size_t i = 0, n = buffer.count(); i < n; ++i)
    {
      require(values[i] >= lowest && values[i] < static_cast<long long>(mapping.size()),
              "invalid source code: " + name);
      if (values[i] >= 0) values[i] = mapping[static_cast<std::size_t>(values[i])];
    }
  }
}

// ---------------------------------------------------------------------------
// Structural check (one chunk's metadata over an entry range)
// ---------------------------------------------------------------------------
struct Check
{
  std::vector<std::string> problems, notes;
  long long entries = 0;
  int groupsChecked = 0;
};

struct IndexRule
{
  const char* name;
  const char* target;  // count branch, or "" for jets / truth jets
  bool allowsNone;
};

const std::vector<IndexRule>& indexRules()
{
  static const std::vector<IndexRule> rules = {
      {"cell_pho", "pho_et", false},  {"isoc_pho", "pho_et", false},
      {"pair_pho", "pho_et", false},  {"pho_assoc_truthpho", "truthpho_pt", true},
      {"plink_reco", "pho_et", true}, {"plink_truth", "truthpho_pt", true},
      {"pair_jet", "", false},        {"jlink_reco", "", true},
      {"jlink_truth", "", true}};
  return rules;
}

json photonIdSpec(const std::string& package, const std::string& system)
{
  const json freeze = readJsonFile(package);
  require(freeze.value("schema", "") == "CanonicalPhotonIDFreezeV1", "not a canonical photon-ID freeze: " + package);
  require(freeze.contains(system), package + " has no " + system + " definition");
  return freeze.at(system).at("features");
}

double bitsToDouble(const json& value)
{
  const auto raw = value.at("f64").get<std::uint64_t>();
  double result;
  std::memcpy(&result, &raw, sizeof(result));
  return result;
}

double valueAt(const Buffer& buffer, std::size_t index)
{
  const unsigned char* p = buffer.bytes() + index * buffer.width;
  if (buffer.type == "Double_t")
  {
    double v;
    std::memcpy(&v, p, 8);
    return v;
  }
  if (buffer.type == "Float_t")
  {
    float v;
    std::memcpy(&v, p, 4);
    return v;
  }
  if (buffer.type == "Int_t")
  {
    std::int32_t v;
    std::memcpy(&v, p, 4);
    return v;
  }
  if (buffer.type == "UInt_t")
  {
    std::uint32_t v;
    std::memcpy(&v, p, 4);
    return v;
  }
  if (buffer.type == "Long64_t")
  {
    std::int64_t v;
    std::memcpy(&v, p, 8);
    return static_cast<double>(v);
  }
  if (buffer.type == "ULong64_t")
  {
    std::uint64_t v;
    std::memcpy(&v, p, 8);
    return static_cast<double>(v);
  }
  throw std::runtime_error("numeric access to " + buffer.name + " (" + buffer.type + ")");
}

long long integerAt(const Buffer& buffer, std::size_t index)
{
  const unsigned char* p = buffer.bytes() + index * buffer.width;
  if (buffer.type == "Int_t")
  {
    std::int32_t v;
    std::memcpy(&v, p, 4);
    return v;
  }
  if (buffer.type == "Long64_t")
  {
    std::int64_t v;
    std::memcpy(&v, p, 8);
    return v;
  }
  if (buffer.type == "UInt_t")
  {
    std::uint32_t v;
    std::memcpy(&v, p, 4);
    return v;
  }
  if (buffer.type == "Short_t")
  {
    std::int16_t v;
    std::memcpy(&v, p, 2);
    return v;
  }
  if (buffer.type == "UShort_t")
  {
    std::uint16_t v;
    std::memcpy(&v, p, 2);
    return v;
  }
  if (buffer.type == "Char_t")
  {
    std::int8_t v;
    std::memcpy(&v, p, 1);
    return v;
  }
  if (buffer.type == "UChar_t")
  {
    std::uint8_t v;
    std::memcpy(&v, p, 1);
    return v;
  }
  if (buffer.type == "ULong64_t")
  {
    std::uint64_t v;
    std::memcpy(&v, p, 8);
    require(v <= static_cast<std::uint64_t>(std::numeric_limits<long long>::max()), buffer.name + " exceeds int64");
    return static_cast<long long>(v);
  }
  throw std::runtime_error("integer access to " + buffer.name + " (" + buffer.type + ")");
}

// The per-chunk structure a reader relies on. `meta` is that chunk's
// metadata (a consolidated file supplies its union codes and counters, and
// its slice defines the entry and retained counts).
Check checkRange(Columns& columns, long long begin, long long end, const json& meta, const eventtree::Schema& schema,
                 const std::string& photonIdPackage, bool consolidatedSlice = false)
{
  Check result;
  const json& info = meta.at("schema");
  const json codes = meta.value("codes", json::object());
  const json collections = meta.value("collections", json{{"reco", json::array()}});
  const long long n = end - begin;
  result.entries = n;
  if (info.value("schema", "") != "RJEventTreeV1") result.problems.push_back("meta/schema is not RJEventTreeV1");
  const std::string written = info.value("schema_sha256", "?");
  if (written != schema.sha256())
    result.notes.push_back("written with schema " + written.substr(0, 12) + "; this checker has " +
                           schema.sha256().substr(0, 12));
  if (!consolidatedSlice && (!info.contains("entries") || info.at("entries") != n))
    result.problems.push_back("entry count " + std::to_string(n) + " differs from meta " +
                              info.value("entries", json()).dump());
  std::vector<std::string> missing;
  for (const auto& branch : schema.branches())
    if (!columns.find(branch.name)) missing.push_back(branch.name);
  if (!missing.empty())
  {
    std::string list;
    for (std::size_t i = 0; i < std::min<std::size_t>(10, missing.size()); ++i) list += (i ? ", " : "") + missing[i];
    result.problems.push_back("branches missing: [" + list + "]");
  }

  const json strides = info.value("strides", json::object());
  std::set<std::string> absent;
  for (const auto& name : info.value("absent_branches", json::array())) absent.insert(name.get<std::string>());

  // Collections in tree order, grouped by schema group (absent excluded).
  std::vector<Buffer*> collectionsInOrder;
  std::map<std::string, std::vector<Buffer*>> groups;
  std::vector<std::string> groupOrder;
  for (const auto& buffer : columns.all())
  {
    if (buffer->name.rfind("n_", 0) == 0) continue;
    if (buffer->kind != Buffer::Kind::Counted && buffer->kind != Buffer::Kind::Strings) continue;
    collectionsInOrder.push_back(buffer.get());
    const auto* branch = schema.find(buffer->name);
    if (!branch || absent.count(buffer->name)) continue;
    if (!groups.count(branch->group)) groupOrder.push_back(branch->group);
    groups[branch->group].push_back(buffer.get());
  }
  struct GroupRule
  {
    Buffer* reference = nullptr;
    std::vector<Buffer*> plain;
    std::vector<std::pair<Buffer*, long long>> strided;
  };
  std::map<std::string, GroupRule> rules;
  for (const auto& group : groupOrder)
  {
    GroupRule rule;
    for (Buffer* buffer : groups[group])
    {
      if (strides.contains(buffer->name))
        rule.strided.emplace_back(buffer, strides.at(buffer->name).get<long long>());
      else if (!rule.reference)
        rule.reference = buffer;
      else
        rule.plain.push_back(buffer);
    }
    if (rule.reference) rules[group] = rule;
  }
  result.groupsChecked = static_cast<int>(rules.size());

  std::vector<std::pair<Buffer*, Buffer*>> counters;
  const json counterMap = info.value("counters", json::object());
  for (const auto& [name, counter] : counterMap.items())
  {
    Buffer* data = columns.find(name);
    Buffer* count = columns.find(counter.get<std::string>());
    if (data && count && (data->kind == Buffer::Kind::Counted || data->kind == Buffer::Kind::Strings))
      counters.emplace_back(data, count);
  }

  std::vector<std::string> recoJets;
  for (const auto& row : collections.value("reco", json::array()))
    recoJets.push_back(row.at("name").get<std::string>() + "_pt");
  struct ActiveIndex
  {
    const IndexRule* rule;
    Buffer* values;
    long long bad = 0;
  };
  std::vector<ActiveIndex> indices;
  for (const auto& rule : indexRules())
  {
    Buffer* values = columns.find(rule.name);
    if (values && values->kind == Buffer::Kind::Counted && !absent.count(rule.name)) indices.push_back({&rule, values});
  }

  struct CodedField
  {
    std::string name;
    Buffer* buffer;
    long long lowest, size;
    bool bad = false;
  };
  std::vector<CodedField> coded;
  for (const auto& [name, dictionary] : codes.items())
    if (Buffer* buffer = columns.find(name))
      coded.push_back({name, buffer, indexField(name) ? -1 : 0, static_cast<long long>(dictionary.size())});

  std::set<std::string> reported;
  auto once = [&](const std::string& message)
  {
    if (reported.insert(message).second) result.problems.push_back(message);
  };
  long long retained = 0, photons = 0;
  bool badRecord = false;
  std::vector<bool> onRejected(collectionsInOrder.size(), false), withValues(collectionsInOrder.size(), false);
  Buffer& record = columns.at("evt_record");
  Buffer* photonEt = columns.find("pho_et");
  for (long long entry = begin; entry < end; ++entry)
  {
    columns.read(entry);
    const long long kind = integerAt(record, 0);
    if (kind != 0 && kind != 1) badRecord = true;
    if (kind == 0) ++retained;
    if (photonEt) photons += static_cast<long long>(photonEt->count());
    for (auto& [group, rule] : rules)
    {
      const std::size_t reference = rule.reference->count();
      for (Buffer* other : rule.plain)
        if (other->count() != reference) once(other->name + " length differs from " + rule.reference->name);
      for (auto& [buffer, stride] : rule.strided)
        if (static_cast<long long>(buffer->count()) != static_cast<long long>(reference) * stride)
          once(buffer->name + " is not " + std::to_string(stride) + " x the " + group + " length");
    }
    for (const auto& [group, per] : std::vector<std::pair<std::string, int>>{{"shv", 7}, {"isow", 2}})
      if (rules.count(group) && rules.count("pho") &&
          rules[group].reference->count() != static_cast<std::size_t>(per) * rules["pho"].reference->count())
        once(group + " is not " + std::to_string(per) + " per photon");
    for (const auto& [data, count] : counters)
      if (static_cast<long long>(data->count()) != integerAt(*count, 0))
        once(count->name + " differs from the length of " + data->name);
    for (auto& index : indices)
    {
      long long limit = 0;
      const std::string name = index.rule->name;
      if (name == "pair_jet" || name == "jlink_reco")
      {
        for (const auto& jets : recoJets)
          if (Buffer* b = columns.find(jets)) limit += static_cast<long long>(b->count());
      }
      else if (name == "jlink_truth")
      {
        for (const char* jets : {"truthjet02_pt", "truthjet03_pt", "truthjet04_pt"})
          if (Buffer* b = columns.find(jets)) limit += static_cast<long long>(b->count());
      }
      else if (Buffer* target = columns.find(index.rule->target))
        limit = static_cast<long long>(target->count());
      const long long lowest = index.rule->allowsNone ? -1 : 0;
      for (std::size_t i = 0, k = index.values->count(); i < k; ++i)
      {
        const long long value = integerAt(*index.values, i);
        if (value < lowest || value >= limit)
        {
          ++index.bad;
          break;
        }
      }
    }
    for (std::size_t i = 0; i < collectionsInOrder.size(); ++i)
      if (collectionsInOrder[i]->count() > 0)
      {
        withValues[i] = true;
        if (kind == 1) onRejected[i] = true;
      }
    for (auto& field : coded)
      for (std::size_t i = 0, k = field.buffer->count(); i < k && !field.bad; ++i)
      {
        const long long value = integerAt(*field.buffer, i);
        if (value < field.lowest || value >= field.size) field.bad = true;
      }
  }
  if (badRecord) result.problems.push_back("evt_record outside {0, 1}");
  if (!consolidatedSlice && (!info.contains("retained_events") || info.at("retained_events") != retained))
    result.problems.push_back("retained entries differ from meta");
  for (const auto& index : indices)
    if (index.bad)
      result.problems.push_back(std::string(index.rule->name) + " out of range in " + std::to_string(index.bad) +
                                " entries");
  // Upstream-rejected entries carry no collection values (first such branch).
  for (std::size_t i = 0; i < collectionsInOrder.size(); ++i)
    if (onRejected[i])
    {
      result.problems.push_back(collectionsInOrder[i]->name + " has values on upstream-rejected entries");
      break;
    }
  for (std::size_t i = 0; i < collectionsInOrder.size(); ++i)
    if (withValues[i] && absent.count(collectionsInOrder[i]->name))
      result.problems.push_back("absent branch " + collectionsInOrder[i]->name + " holds values");
  for (const auto& field : coded)
    if (field.bad) result.problems.push_back(field.name + " code outside its dictionary");
  if (!result.problems.empty())
  {
    result.notes.push_back("photon-ID inputs not checked: the file has structure problems");
    return result;
  }
  const json constants = meta.value("constants", json::object());
  if (photonIdPackage.empty() || !photons || !constants.contains("showerViews.definition_name")) return result;

  // Photon-ID model inputs: one slot per shower view, the model's view with
  // the frozen semantics and floor, and stored features equal to the named
  // branches (compared as float32, NaN equal to NaN).
  const std::string system = info.value("system", "");
  const json spec = photonIdSpec(photonIdPackage, system);
  const auto names = constants.at("showerViews.definition_name").get<std::vector<std::string>>();
  const std::string view = spec.at("view").get<std::string>();
  const auto occurrences = std::count(names.begin(), names.end(), view);
  if (occurrences != 1)
  {
    result.problems.push_back("shower view " + view + " appears " + std::to_string(occurrences) +
                              " times in the view order");
    return result;
  }
  require(names.size() == 7, "the shower-view order must hold the seven stored views");
  const std::size_t position = static_cast<std::size_t>(std::find(names.begin(), names.end(), view) - names.begin());
  if (constants.at("showerViews.semantic_sha256").at(position) != spec.at("view_semantic_sha256"))
    result.problems.push_back(view + " view semantic hash differs from the frozen model contract");
  if (bitsToDouble(constants.at("showerViews.floor_gev").at(position)) != spec.at("floor_gev").get<double>())
    result.problems.push_back(view + " view floor differs from the frozen model contract");
  const long long stride = info.at("strides").at("shv_features").get<long long>();
  const auto order = spec.at("order").get<std::vector<std::string>>();
  if (stride != spec.at("count").get<long long>())
  {
    result.problems.push_back(view + " view carries " + std::to_string(stride) + " features, the model " +
                              std::to_string(spec.at("count").get<long long>()));
    return result;
  }
  std::vector<std::string> inputs;
  for (const auto& name : order) inputs.push_back(eventtree::featureBranch(name));
  long long mismatches = 0;
  bool fractionalCentrality = false;
  for (long long entry = begin; entry < end; ++entry)
  {
    columns.read(entry);
    const std::size_t count = photonEt->count();
    const Buffer& features = columns.at("shv_features");
    for (std::size_t photon = 0; photon < count; ++photon)
    {
      bool same = true;
      for (std::size_t input = 0; input < inputs.size(); ++input)
      {
        const std::string& branch = inputs[input];
        const Buffer& source = columns.at(branch);
        double value;
        if (branch.rfind("shv_", 0) == 0)
          value = valueAt(source, photon * 7 + position);
        else if (branch == "vtx_z" || branch == "cent_value")
          value = valueAt(source, 0);
        else
          value = valueAt(source, photon);
        const float rebuilt = static_cast<float>(value);
        const float stored =
            static_cast<float>(valueAt(features, (photon * 7 + position) * static_cast<std::size_t>(stride) + input));
        if (!(stored == rebuilt || (std::isnan(stored) && std::isnan(rebuilt)))) same = false;
      }
      mismatches += same ? 0 : 1;
      if (system == "auau")
      {
        const double centrality = valueAt(columns.at("cent_value"), 0);
        if (std::isfinite(centrality) && centrality != std::nearbyint(centrality)) fractionalCentrality = true;
      }
    }
  }
  if (mismatches)
    result.problems.push_back(std::to_string(mismatches) + " photons: " + view +
                              " features differ from the named shower branches");
  if (fractionalCentrality) result.problems.push_back("Au+Au centrality is not an integer bin for every photon");
  result.notes.push_back("photon-ID inputs checked on " + std::to_string(photons) + " photons (" + view + " view)");
  return result;
}

// ---------------------------------------------------------------------------
// Source assembly helpers
// ---------------------------------------------------------------------------
void expectedRanges(const json& chunks, long long expectedEnd, const json& exclusions)
{
  std::vector<std::pair<long long, long long>> ranges;
  for (const auto& list : {chunks, exclusions})
    for (const auto& item : list)
      ranges.emplace_back(item.at("event_begin").get<long long>(), item.at("event_end").get<long long>());
  std::sort(ranges.begin(), ranges.end());
  long long end = 0;
  for (const auto& [begin, stop] : ranges)
  {
    require(begin == end && stop >= end, "source range gap, overlap or ordering error");
    end = stop;
  }
  require(end == expectedEnd, "source range does not reach expected end");
}

// Union of observed shapes. Per-chunk metadata remains the inverse-proof
// authority; zero strides carry no observation, conflicts refuse.
std::pair<json, json> mergeShapes(const std::vector<json>& metas)
{
  json strides = json::object();
  json collections = metas.front().at("collections");
  for (const json& meta : metas)
  {
    for (const auto& [name, value] : meta.at("schema").at("strides").items())
    {
      require(value.is_number_integer() && value.get<long long>() >= 0, "invalid stride: " + name);
      const long long old = strides.value(name, 0LL), now = value.get<long long>();
      require(!old || !now || old == now, "conflicting nonzero stride: " + name);
      strides[name] = old ? old : now;
    }
    std::set<std::string> a, b;
    for (const auto& item : meta.at("collections").items()) a.insert(item.key());
    for (const auto& item : collections.items()) b.insert(item.key());
    require(a == b, "jet collection families differ");
    for (const auto& [family, rows] : meta.at("collections").items())
    {
      json& expected = collections[family];
      require(rows.size() == expected.size(), "jet collection roster differs");
      for (std::size_t i = 0; i < rows.size(); ++i)
      {
        json& known = expected[i];
        const json& row = rows[i];
        for (const char* key : {"name", "input_identity", "radius", "radius_code"})
          require(known.value(key, json()) == row.value(key, json()),
                  std::string("jet collection identity differs: ") + key);
        for (const auto& [key, value] : row.items())
        {
          require(!known.contains(key) || known.at(key) == value, "jet collection recipe differs: " + key);
          known[key] = value;
        }
      }
    }
  }
  for (const json& meta : metas)
    for (const auto& item : strides.items())
    {
      const json& absent = meta.at("schema").at("absent_branches");
      require(meta.at("schema").at("strides").contains(item.key()) ||
                  std::find(absent.begin(), absent.end(), item.key()) != absent.end(),
              "non-absent branch lacks stride: " + item.key());
    }
  return {strides, collections};
}

// Physics identity of a chunk, independent of its statistics and ROOT build.
json definition(const json& meta, const eventtree::Schema& schema)
{
  require(schema.acceptsDigest(meta.at("schema").at("schema_sha256").get<std::string>()), "source schema differs");
  const json id = meta.value("photon_id", json::object()), iso = meta.value("isolation", json::object());
  return json::array({schema.contractSha256(), id.value("freeze_sha256", json()), id.value("model_sha256", json()),
                      iso.value("freeze_sha256", json())});
}

std::map<std::string, std::vector<int>> codeMaps(const json& chunk)
{
  std::map<std::string, std::vector<int>> maps;
  const json mappings = chunk.value("code_maps", json::object());
  for (const auto& [name, mapping] : mappings.items()) maps[name] = mapping.get<std::vector<int>>();
  return maps;
}

// Carried objects: every non-tree, non-directory key (highest cycle).
std::vector<std::pair<std::string, TKey*>> carriedObjects(TFile& file, bool directMetadata)
{
  std::vector<std::pair<std::string, TKey*>> objects;
  visitKeys(&file, "",
            [&](const std::string& path, TKey* key)
            {
              if (isTree(key)) return;
              const bool direct = path.rfind("meta/", 0) == 0 && path.find('/', 5) == std::string::npos &&
                                  std::string(key->GetClassName()) == "TObjString";
              if (direct && !directMetadata) return;
              objects.emplace_back(path, key);
            });
  return objects;
}

std::string checkReportJson(const Check& check, const std::string& path)
{
  return json{{"ok", check.problems.empty()},
              {"entries", check.entries},
              {"problems", check.problems},
              {"notes", check.notes},
              {"groups_checked", check.groupsChecked},
              {"path", path}}
      .dump();
}

Check checkConsolidated(TFile& file, const eventtree::Schema& schema, const std::string& photonIdPackage)
{
  const json meta = readMetadata(file);
  TTree* tree = eventTree(file);
  const json& chunks = meta.at("consolidation").at("chunks");
  const json& consolidation = meta.at("consolidation");
  expectedRanges(chunks,
                 consolidation.value("original_source_events", consolidation.at("expected_events").get<long long>()),
                 consolidation.value("exclusions", json::array()));
  std::vector<json> originals;
  for (const auto& chunk : chunks) originals.push_back(chunk.at("original_meta"));
  const auto [strides, collections] = mergeShapes(originals);
  require(strides == meta.at("schema").at("strides") && collections == meta.at("collections"),
          "consolidated shape union differs");
  long long end = 0;
  for (const auto& chunk : chunks)
  {
    require(chunk.at("entry_begin") == end && chunk.at("entry_end").get<long long>() >= end,
            "invalid consolidated chunk coverage");
    require(chunk.at("entry_end").get<long long>() - end ==
                chunk.at("original_meta").at("schema").at("entries").get<long long>(),
            "chunk accounting metadata differs");
    end = chunk.at("entry_end").get<long long>();
  }
  require(end == tree->GetEntries() && end == consolidation.at("expected_events").get<long long>(),
          "consolidated manifest coverage differs");
  Columns columns(tree);
  Check total;
  for (const auto& chunk : chunks)
  {
    json local = chunk.at("original_meta");
    local["codes"] = meta.at("codes");
    local["schema"]["counters"] = meta.at("schema").at("counters");
    const long long begin = chunk.at("entry_begin").get<long long>(), stop = chunk.at("entry_end").get<long long>();
    const Check part = checkRange(columns, begin, stop, local, schema, photonIdPackage, true);
    total.problems.insert(total.problems.end(), part.problems.begin(), part.problems.end());
    total.notes.insert(total.notes.end(), part.notes.begin(), part.notes.end());
    total.entries += part.entries;
    total.groupsChecked = std::max(total.groupsChecked, part.groupsChecked);
  }
  if (total.entries != tree->GetEntries()) total.problems.push_back("chunk entry coverage differs");
  if (treeKeys(file).size() != 1) total.problems.push_back("not exactly one TTree");
  return total;
}

json proveConsolidatedFile(const std::string& path, const eventtree::Schema& schema, const std::string& photonIdPackage)
{
  long long checked = 0, objects = 0;
  auto final = openRead(path);
  const json meta = readMetadata(*final);
  const json& chunks = meta.at("consolidation").at("chunks");
  TTree* tree = eventTree(*final);
  require(tree->GetEntries() == meta.at("consolidation").at("expected_events").get<long long>(),
          "consolidated entry count differs");
  Columns output(tree);
  for (const auto& chunk : chunks)
  {
    const std::string piece = chunk.at("path").get<std::string>();
    require(eventtree::fileSha256(piece) == chunk.at("sha256").get<std::string>(), "source hash drift during proof");
    auto source = openRead(piece);
    require(readMetadata(*source) == chunk.at("original_meta"), "chunk provenance changed");
    Columns input(eventTree(*source));
    const auto maps = codeMaps(chunk);
    const long long offset = chunk.at("entry_begin").get<long long>();
    for (long long entry = 0; entry < input.entries(); ++entry)
    {
      input.read(entry);
      remap(input, maps);
      output.read(offset + entry);
      for (const auto& buffer : input.all())
      {
        const Buffer* saved = output.find(buffer->name);
        require(saved && sameValues(*buffer, *saved), "consolidated payload differs: " + buffer->name);
      }
      ++checked;
    }
    const std::string prefix = chunk.at("object_prefix").get<std::string>();
    for (const auto& [name, key] : carriedObjects(*source, true))
    {
      TKey* copy = latestKeyAt(*final, prefix + name);
      require(copy != nullptr, "chunk object missing: " + name);
      auto before = readObject(key);
      auto after = readObject(copy);
      require(streamed(*before) == streamed(*after), "chunk object differs: " + name);
      ++objects;
    }
    require(eventtree::fileSha256(piece) == chunk.at("sha256").get<std::string>(), "input changed during proof");
  }
  require(checked == tree->GetEntries(), "conservation coverage differs");
  const Check structure = checkConsolidated(*final, schema, photonIdPackage);
  require(structure.problems.empty(),
          "consolidated structure: " + (structure.problems.empty() ? std::string() : structure.problems.front()));
  return {{"status", "PASS"}, {"entries", checked}, {"objects", objects}, {"chunks", chunks.size()}};
}

// ---------------------------------------------------------------------------
// Group publication helpers
// ---------------------------------------------------------------------------
fs::path safeRelative(const std::string& name)
{
  const fs::path path(name);
  bool parent = false;
  for (const auto& part : path) parent |= part == "..";
  require(!path.is_absolute() && !parent && path.extension() == ".root", "unsafe logical output name: " + name);
  return path;
}

// Unique name for a temporary or candidate file in a shared directory.
std::string newToken()
{
  std::ifstream random("/dev/urandom", std::ios::binary);
  std::array<unsigned char, 16> bytes{};
  random.read(reinterpret_cast<char*>(bytes.data()), bytes.size());
  require(random.gcount() == static_cast<std::streamsize>(bytes.size()), "cannot read /dev/urandom");
  std::ostringstream out;
  for (const unsigned char byte : bytes) out << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(byte);
  return out.str();
}

// DATA physical events must be unique across a group. Keys are sorted in
// bounded runs on disk and merged, so memory does not grow with the group.
json populationCheck(const json& files, const fs::path& scratch)
{
  constexpr std::size_t kRun = 1 << 22;  // keys per sorted run (64 MB)
  using Key = std::pair<long long, long long>;
  json reference;
  long long entries = 0, physical = 0;
  std::vector<fs::path> runs;
  std::vector<Key> block;
  auto flush = [&]()
  {
    if (block.empty()) return;
    std::sort(block.begin(), block.end());
    const fs::path run = scratch / ("run-" + std::to_string(runs.size()) + ".bin");
    std::ofstream out(run, std::ios::binary | std::ios::trunc);
    out.write(reinterpret_cast<const char*>(block.data()), static_cast<std::streamsize>(block.size() * sizeof(Key)));
    require(static_cast<bool>(out), "cannot write population index " + run.string());
    runs.push_back(run);
    block.clear();
  };
  fs::create_directories(scratch);
  for (const auto& item : files)
  {
    auto file = openRead(item.at("path").get<std::string>());
    const json source = readMetadata(*file).at("source");
    const std::string role = source.value("si_di_role", "");
    const bool data = role == "DATA";
    const json identity = json::array({source.value("dataset", json()), source.value(data ? "run" : "sample", json()),
                                       source.value("period", json()), source.value("si_di_role", json())});
    if (reference.is_null()) reference = identity;
    require(identity == reference, "group mixes DATA runs or SIM sample identity");
    require(source.value("logical_source_id", json()) == item.at("logical_id"),
            "file logical source differs from group manifest");
    TTree* tree = eventTree(*file);
    Int_t run = 0, valid = 0;
    Long64_t sequence = 0;
    tree->SetBranchStatus("*", false);
    for (const char* name : {"evt_run", "evt_phys_seq", "evt_phys_seq_valid"}) tree->SetBranchStatus(name, true);
    tree->SetBranchAddress("evt_run", &run);
    tree->SetBranchAddress("evt_phys_seq", &sequence);
    tree->SetBranchAddress("evt_phys_seq_valid", &valid);
    for (Long64_t entry = 0, n = tree->GetEntries(); entry < n; ++entry)
    {
      require(tree->GetEntry(entry) > 0, "cannot read event identity");
      ++entries;
      if (!data) continue;
      require(valid == 1, "DATA group has unverifiable physical event identity");
      block.emplace_back(run, sequence);
      ++physical;
      if (block.size() == kRun) flush();
    }
  }
  flush();
  struct Stream
  {
    std::ifstream in;
    Key key;
  };
  std::vector<std::unique_ptr<Stream>> streams;
  auto next = [](Stream& s) { return static_cast<bool>(s.in.read(reinterpret_cast<char*>(&s.key), sizeof(Key))); };
  auto later = [](const Stream* a, const Stream* b) { return a->key > b->key; };
  std::priority_queue<Stream*, std::vector<Stream*>, decltype(later)> heap(later);
  for (const auto& run : runs)
  {
    streams.push_back(std::make_unique<Stream>());
    streams.back()->in.open(run, std::ios::binary);
    if (next(*streams.back())) heap.push(streams.back().get());
  }
  bool first = true;
  Key previous{};
  while (!heap.empty())
  {
    Stream* s = heap.top();
    heap.pop();
    require(first || s->key != previous, "duplicate physical DATA event across logical sources");
    previous = s->key;
    first = false;
    if (next(*s)) heap.push(s);
  }
  streams.clear();
  for (const auto& run : runs) fs::remove(run);
  return {
      {"status", "PASS"}, {"entries", entries}, {"unique_data_events", physical}, {"logical_sources", files.size()}};
}

json resolveGroupJson(const std::string& group)
{
  const fs::path root = fs::absolute(group).lexically_normal();
  const fs::path current = root / "CURRENT";
  require(fs::is_symlink(current), "group is not published: " + root.string());
  const fs::path version = fs::canonical(current);
  require(version.parent_path() == fs::canonical(root) / "versions",
          "publication points outside immutable version root");
  const json manifest = readJsonFile(version / "MANIFEST.json");
  json body = manifest;
  body.erase("version");
  require(eventtree::textSha256(dumps(body)) == manifest.at("version").get<std::string>() &&
              manifest.at("version").get<std::string>() == version.filename().string(),
          "group manifest changed");
  return {{"directory", version.string()}, {"manifest", manifest}};
}
}  // namespace

// ===========================================================================
// Public operations
// ===========================================================================
namespace eventtree
{
std::string metadataJson(const std::string& path)
{
  auto file = openRead(path);
  return readMetadata(*file).dump();
}

std::string checkFile(const std::string& path, const std::string& schemaPath, const std::string& photonIdPackage)
{
  const Schema& schema = Schema::load(schemaPath);
  auto file = openRead(path);
  const json meta = readMetadata(*file);
  if (meta.contains("consolidation")) return checkReportJson(checkConsolidated(*file, schema, photonIdPackage), path);
  Columns columns(eventTree(*file));
  Check check = checkRange(columns, 0, columns.entries(), meta, schema, photonIdPackage);
  const auto trees = treeKeys(*file);
  if (trees.size() != 1 || trees.front().rfind(std::string(kTreeName) + ";", 0) != 0)
  {
    std::string list;
    for (const auto& key : trees) list += (list.empty() ? "" : ", ") + key;
    check.problems.push_back("the file holds " + std::to_string(trees.size()) + " TTree keys, not one eventTree: [" +
                             list + "]");
  }
  return checkReportJson(check, path);
}

std::string checkProducerOutput(const std::string& path, const std::string& schemaPath, const std::string& jobJson)
{
  const Schema& schema = Schema::load(schemaPath);
  const json job = parse(jobJson, "producer job");
  auto file = openRead(path);
  const auto completion = keyValues(*file, "completion");
  const auto execution = keyValues(*file, "metadata");
  require(completion.count("completion_status") && completion.at("completion_status") == "complete",
          "producer output is not complete");
  auto text = [](const json& value) { return value.is_string() ? value.get<std::string>() : value.dump(); };
  std::map<std::string, std::string> expected;
  for (const char* key : {"run", "segment", "source_file_ordinal", "first_entry", "source_manifest_sha256", "sample",
                          "period", "si_di_role", "configuration_sha256", "macro_sha256"})
    expected[key] = text(job.at(key));
  expected["max_events"] = text(job.at("events"));
  expected["producer_library_sha256"] = job.at("producer_library_sha256").get<std::string>();
  for (const auto& [key, value] : expected)
    require(execution.count(key) && execution.at(key) == value, "output provenance mismatch: " + key);
  const json meta = readMetadata(*file);
  require(meta.at("schema").at("schema_sha256") == schema.sha256(), "native schema differs");
  require(meta.at("native_capture").at("missing_branches").empty(), "native mapping is incomplete");
  require(meta.at("source").value("completed", false) == true, "source is incomplete");
  for (const char* key : {"run", "segment", "source_file_ordinal", "first_entry", "sample", "period", "si_di_role",
                          "source_manifest_sha256"})
    require(meta.at("source").value(key, json()) == job.at(key), std::string("source metadata mismatch: ") + key);
  const auto trees = treeKeys(*file);
  require(trees.size() == 1 && trees.front().rfind(std::string(kTreeName) + ";", 0) == 0,
          "producer must write exactly one eventTree");
  TTree* tree = eventTree(*file);
  for (const auto& branch : schema.branches())
    require(tree->GetBranch(branch.name.c_str()) != nullptr, "producer omitted schema branch " + branch.name);
  const long long observed = std::stoll(completion.at("upstream_observed_events"));
  const long long retained = std::stoll(completion.at("retained_events"));
  const long long rejected = std::stoll(completion.at("upstream_rejected_events"));
  require(observed == job.at("events").get<long long>(), "source ended before/after its declared event count");
  require(std::stoll(completion.at("encountered_events")) == retained && observed == retained + rejected,
          "producer completion accounting is inconsistent");
  require(tree->GetEntries() == observed && meta.at("schema").at("entries") == observed, "saved event count differs");
  require(meta.at("schema").at("retained_events") == retained, "retained metadata count differs");
  Int_t record = 0, run = 0;
  Long64_t ordinal = 0, entry = 0;
  tree->SetBranchStatus("*", false);
  for (const char* name : {"evt_record", "evt_run", "evt_src_file_ordinal", "evt_src_entry"})
    tree->SetBranchStatus(name, true);
  tree->SetBranchAddress("evt_record", &record);
  tree->SetBranchAddress("evt_run", &run);
  tree->SetBranchAddress("evt_src_file_ordinal", &ordinal);
  tree->SetBranchAddress("evt_src_entry", &entry);
  long long kept = 0, dropped = 0;
  const long long first = job.at("first_entry").get<long long>();
  const Long64_t entries = tree->GetEntries();
  for (Long64_t i = 0; i < entries; ++i)
  {
    require(tree->GetEntry(i) > 0, "cannot read saved event accounting");
    require(entry == first + i, "source entries contain a gap, overlap or reordered event");
    require(ordinal == job.at("source_file_ordinal").get<long long>(), "wrong source ordinal");
    require(run == job.at("run").get<long long>(), "wrong event run");
    require(record == 0 || record == 1, "unknown event record kind");
    (record == 0 ? kept : dropped) += 1;
  }
  require(kept == retained && dropped == rejected, "saved event accounting differs");
  json result = json::object();
  for (const auto& [key, value] : completion) result[key] = value;
  return result.dump();
}

std::string consolidate(const std::string& requestJson)
{
  const json request = parse(requestJson, "consolidation request");
  const Schema& schema = Schema::load(request.at("schema_path").get<std::string>());
  const std::string package = request.value("photon_id_package", "");
  const json& chunks = request.at("chunks");
  const std::string output = request.at("output").get<std::string>();
  const std::string identity = request.at("logical_source_id").get<std::string>();
  const long long expectedEvents = request.at("expected_events").get<long long>();
  const json exclusions = request.value("exclusions", json::array());
  require(!chunks.empty(), "no source chunks");
  std::set<std::string> hashes, pieces;
  for (const auto& chunk : chunks)
  {
    require(hashes.insert(chunk.at("sha256").get<std::string>()).second, "same payload assigned to multiple chunks");
    require(pieces.insert(chunk.at("piece_id").get<std::string>()).second, "duplicate piece identity");
  }
  for (const auto& exclusion : exclusions)
    require(!exclusion.value("authority_json", std::string()).empty(), "excluded ranges lack pinned authority");
  expectedRanges(chunks, expectedEvents, exclusions);
  long long admitted = 0;
  for (const auto& chunk : chunks)
    admitted += chunk.at("event_end").get<long long>() - chunk.at("event_begin").get<long long>();

  // Read and bind every piece before creating the candidate.
  std::vector<json> metas;
  json unions = json::object();
  long long entries = 0, retained = 0;
  std::vector<std::pair<bool, long long>> centrality;
  for (const auto& chunk : chunks)
  {
    const std::string path = chunk.at("path").get<std::string>();
    require(fileSha256(path) == chunk.at("sha256").get<std::string>(), "input changed: " + path);
    auto file = openRead(path);
    json meta = readMetadata(*file);
    TTree* tree = eventTree(*file);
    centrality.emplace_back(latestKeyAt(*file, "hist/SIM/h_centrality") != nullptr, tree->GetEntries());
    require(!meta.contains("consolidation"), "source is already consolidated");
    require(schema.acceptsDigest(meta.at("schema").at("schema_sha256").get<std::string>()), "source schema differs");
    require(tree->GetEntries() == meta.at("schema").at("entries").get<long long>(), "source entry metadata differs");
    require(meta.at("schema").at("entries").get<long long>() ==
                chunk.at("event_end").get<long long>() - chunk.at("event_begin").get<long long>(),
            "source event accounting incomplete: retained plus rejected must cover interval");
    if (!metas.empty())
    {
      require(definition(meta, schema) == definition(metas.front(), schema), "nominal definition differs");
      for (const char* key : {"system", "lane"})
        require(meta.at("schema").value(key, json()) == metas.front().at("schema").value(key, json()),
                std::string("chunk schema property differs: ") + key);
      for (const char* key : {"dataset", "run", "segment", "sample", "period", "si_di_role"})
        require(meta.at("source").value(key, json()) == metas.front().at("source").value(key, json()),
                std::string("chunk source differs: ") + key);
    }
    const json pieceCodes = meta.value("codes", json::object());
    for (const auto& [branch, labels] : pieceCodes.items())
    {
      // Each chunk code must map to exactly one union entry and back.
      std::set<std::string> distinct;
      for (const auto& label : labels) distinct.insert(label.dump());
      require(distinct.size() == labels.size(), "ambiguous code dictionary in " + path + ": " + branch);
      json& target = unions[branch];
      if (target.is_null()) target = json::array();
      for (const auto& label : labels)
        if (std::find(target.begin(), target.end(), label) == target.end()) target.push_back(label);
    }
    entries += tree->GetEntries();
    retained += meta.at("schema").at("retained_events").get<long long>();
    metas.push_back(std::move(meta));
  }
  require(entries == admitted, "event count differs from admitted source coverage");
  if (std::any_of(centrality.begin(), centrality.end(), [](const auto& c) { return c.first; }))
    require(std::all_of(centrality.begin(), centrality.end(), [](const auto& c) { return c.first || c.second == 0; }),
            "nonempty chunk is missing the source centrality histogram");

  json combined = metas.front();
  combined["codes"] = unions;
  const auto [strides, collections] = mergeShapes(metas);
  combined["schema"]["strides"] = strides;
  combined["collections"] = collections;
  std::set<std::string> absent;
  for (const auto& name : metas.front().at("schema").at("absent_branches")) absent.insert(name.get<std::string>());
  for (const json& meta : metas)
  {
    std::set<std::string> here;
    for (const auto& name : meta.at("schema").at("absent_branches")) here.insert(name.get<std::string>());
    std::set<std::string> both;
    std::set_intersection(absent.begin(), absent.end(), here.begin(), here.end(), std::inserter(both, both.begin()));
    absent = both;
  }
  combined["schema"]["schema_sha256"] = schema.sha256();
  combined["schema"]["entries"] = entries;
  combined["schema"]["retained_events"] = retained;
  combined["schema"]["absent_branches"] = absent;
  // Cumulative run trigger tables are kept per chunk, never summed.
  const json trigger = metas.front().value("trigger_run_info", json());
  if (!std::all_of(metas.begin(), metas.end(),
                   [&](const json& m) { return m.value("trigger_run_info", json()) == trigger; }))
    combined.erase("trigger_run_info");
  if (trigger.is_object() && trigger.contains("name"))
  {
    const json selection = {{"name", trigger.at("name")}, {"initial_prescale", trigger.at("initial_prescale")}};
    for (const json& meta : metas)
    {
      const json other = meta.value("trigger_run_info", json());
      require(other.is_object() && other.contains("name") &&
                  json({{"name", other.at("name")}, {"initial_prescale", other.value("initial_prescale", json())}}) ==
                      selection,
              "trigger names/prescales vary within logical source");
    }
    combined["trigger_selection"] = selection;
  }
  combined["source"]["logical_source_id"] = identity;
  combined["source"]["consolidated"] = true;
  combined["consolidation"] = {
      {"schema", "RJLogicalSourceV1"},
      {"logical_source_id", identity},
      {"expected_events", admitted},
      {"original_source_events", expectedEvents},
      {"exclusions", exclusions},
      {"chunks", json::array()},
      {"metadata_policy", "per-chunk originals; cumulative counters and normalization are not added"},
      {"content_proof", "payload bitwise after reversible dictionary mapping"}};

  Progress progress(request.value("progress_path", ""), "source_validation", identity, entries);
  {
    NoHistogramDirectory guard;
    auto out = createFile(output);
    TTree* tree = nullptr;
    long long offset = 0;
    std::vector<std::pair<long long, long long>> physical;
    std::map<std::string, std::string> layout;
    json counters = json::object();
    for (std::size_t index = 0; index < chunks.size(); ++index)
    {
      const json& piece = chunks[index];
      const json& meta = metas[index];
      std::map<std::string, std::vector<int>> maps;
      json mapJson = json::object();
      const json pieceCodes = meta.value("codes", json::object());
      for (const auto& [branch, labels] : pieceCodes.items())
      {
        std::vector<int> mapping;
        const json& target = unions.at(branch);
        for (const auto& label : labels)
          mapping.push_back(static_cast<int>(std::find(target.begin(), target.end(), label) - target.begin()));
        maps[branch] = mapping;
        mapJson[branch] = mapping;
      }
      std::ostringstream prefix;
      prefix << "meta/chunks/" << std::setw(6) << std::setfill('0') << index << "/objects/";
      json record = piece;
      record["entry_begin"] = offset;
      record["entry_end"] = offset + meta.at("schema").at("entries").get<long long>();
      record["original_meta"] = meta;
      record["code_maps"] = mapJson;
      record["object_prefix"] = prefix.str();
      combined["consolidation"]["chunks"].push_back(record);

      auto source = openRead(piece.at("path").get<std::string>());
      TTree* input = eventTree(*source);
      Columns columns(input);
      out->cd();
      if (layout.empty()) layout = columns.layout();
      require(columns.layout() == layout, "chunk branch layout differs: " + piece.at("path").get<std::string>());
      if (!tree)
      {
        tree = input->CloneTree(0);
        require(tree != nullptr, "cannot book the consolidated eventTree");
        tree->SetDirectory(out.get());
        tree->SetAutoFlush(-16 * 1024 * 1024);
      }
      columns.attach(tree);
      Buffer& sourceEntry = columns.at("evt_src_entry");
      Buffer& run = columns.at("evt_run");
      Buffer& sequence = columns.at("evt_phys_seq");
      Buffer& valid = columns.at("evt_phys_seq_valid");
      const long long begin = piece.at("event_begin").get<long long>();
      for (long long entry = 0; entry < columns.entries(); ++entry)
      {
        columns.read(entry);
        require(integerAt(sourceEntry, 0) == begin + entry, "missing, duplicated or wrongly associated source events");
        if (integerAt(valid, 0) == 1) physical.emplace_back(integerAt(run, 0), integerAt(sequence, 0));
        remap(columns, maps);
        out->cd();
        require(tree->Fill() >= 0, "TTree::Fill failed");
        progress.update(offset + entry + 1);
      }
      offset = record.at("entry_end").get<long long>();
      tree->ResetBranchAddresses();
      require(fileSha256(piece.at("path").get<std::string>()) == piece.at("sha256").get<std::string>(),
              "input changed during consolidation");
    }
    std::sort(physical.begin(), physical.end());
    require(std::adjacent_find(physical.begin(), physical.end()) == physical.end(),
            "duplicate physical event within logical source");
    for (TObject* item : *tree->GetListOfBranches())
    {
      auto* branch = static_cast<TBranch*>(item);
      if (branch->GetListOfLeaves()->GetEntries() != 1) continue;
      if (TLeaf* count = static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0))->GetLeafCount())
        counters[branch->GetName()] = count->GetName();
    }
    combined["schema"]["counters"] = counters;
    combined["schema"]["encoding"] = "counted";
    out->cd();
    require(tree->Write("", TObject::kOverwrite) > 0, "TTree::Write failed");
    for (const auto& [key, value] : combined.items()) writeText(*out, "meta/" + key, dumps(value));
    // Per-chunk originals: every carried object, byte preserved, under its prefix.
    std::unique_ptr<TH1> aggregate;
    for (std::size_t index = 0; index < chunks.size(); ++index)
    {
      auto source = openRead(chunks[index].at("path").get<std::string>());
      const std::string prefix = combined["consolidation"]["chunks"][index]["object_prefix"].get<std::string>();
      for (const auto& [name, key] : carriedObjects(*source, true))
      {
        auto object = readObject(key);
        writeObject(*out, prefix + name, *object);
        // The additive SIM centrality denominator for the maintained reader.
        if (name == "hist/SIM/h_centrality")
        {
          auto* histogram = dynamic_cast<TH1*>(object.get());
          require(histogram != nullptr, "centrality object is not a histogram");
          if (!aggregate)
          {
            aggregate.reset(static_cast<TH1*>(histogram->Clone("h_centrality")));
            aggregate->SetDirectory(nullptr);
          }
          else
          {
            bool edges = aggregate->GetNbinsX() == histogram->GetNbinsX();
            for (int bin = 1; edges && bin <= histogram->GetNbinsX() + 1; ++bin)
              edges = aggregate->GetXaxis()->GetBinLowEdge(bin) == histogram->GetXaxis()->GetBinLowEdge(bin);
            require(edges, "centrality histogram edges differ");
            require(aggregate->Add(histogram), "centrality histogram merge failed");
          }
        }
      }
    }
    if (aggregate) writeObject(*out, "hist/SIM/h_centrality", *aggregate);
    require(!out->TestBit(TFile::kWriteError), "ROOT file write error");
    out->Close();
  }
  const json proof = proveConsolidatedFile(output, schema, package);
  return json{
      {"schema", "RJLogicalSourceReceiptV1"}, {"status", "PASS"},   {"logical_source_id", identity}, {"output", output},
      {"output_sha256", fileSha256(output)},  {"entries", entries}, {"pieces", chunks.size()},       {"proof", proof}}
      .dump();
}

std::string proveConsolidated(const std::string& path, const std::string& schemaPath)
{
  return proveConsolidatedFile(path, Schema::load(schemaPath), "").dump();
}

// ===========================================================================
// Branch additions
// ===========================================================================
struct BranchAddition::State
{
  json request;
  std::string original, candidate;
  std::unique_ptr<TFile> source, output;
  TTree* input = nullptr;
  TTree* tree = nullptr;
  TLeaf* photonCount = nullptr;
  struct Added
  {
    std::string name, type;
    Int_t count = 0;
    std::vector<double> values;
    std::vector<long long> counts;
    std::size_t consumed = 0, valueOffset = 0;
    std::vector<double> buffer;
  };
  std::vector<Added> added;
  long long next = 0;
  bool finished = false;
};

BranchAddition::BranchAddition() = default;

void BranchAddition::open(const std::string& requestJson)
{
  require(!m_state, "branch addition already opened");
  auto state = std::make_unique<State>();
  State& s = *state;
  s.request = parse(requestJson, "branch addition request");
  s.original = s.request.at("original").get<std::string>();
  s.candidate = s.request.at("candidate").get<std::string>();
  require(fs::absolute(s.original).lexically_normal() != fs::absolute(s.candidate).lexically_normal(),
          "update must use a separate candidate");
  require(fileSha256(s.original) == s.request.at("input_sha256").get<std::string>(), "update input hash changed");
  const Schema& schema = Schema::load(s.request.at("schema_path").get<std::string>());
  s.source = openRead(s.original);
  s.input = eventTree(*s.source);
  const json meta = readMetadata(*s.source);
  std::set<std::string> payload, declared;
  for (TObject* item : *s.input->GetListOfBranches())
    if (std::string name = item->GetName(); name.rfind("n_", 0) != 0) payload.insert(name);
  for (const auto& branch : schema.branches()) declared.insert(branch.name);
  const json models = meta.value("photon_models", json::object());
  for (const auto& [name, receipt] : models.items()) declared.insert(receipt.at("score_branch").get<std::string>());
  for (const auto& item : s.request.at("added"))
  {
    State::Added added;
    added.name = item.at("name").get<std::string>();
    added.type = item.value("type", "double");
    require(added.type == "double", "added branches are vector<double> collections: " + added.name);
    require(added.name.rfind("pho_", 0) == 0 && !schema.find(added.name),
            "optional branch collides with base schema: " + added.name);
    require(!payload.count(added.name), "addition would overwrite an existing branch: " + added.name);
    require(declared.insert(added.name).second, "score branch already exists: " + added.name);
    s.added.push_back(std::move(added));
  }
  std::set<std::string> after = payload;
  for (const auto& added : s.added) after.insert(added.name);
  require(after == declared, "update schema has missing or undeclared payloads");
  // Keep every old basket and metadata key at its original byte offset. The
  // candidate is a create-once copy; UPDATE appends only the new score data.
  // Re-filling old branches would change basket packing and duplicate work.
  const fs::path candidate(s.candidate);
  if (candidate.has_parent_path()) fs::create_directories(candidate.parent_path());
  require(fs::copy_file(s.original, candidate, fs::copy_options::none), "cannot copy augmentation input");
  fs::permissions(candidate, fs::perms::owner_write, fs::perm_options::add);
  s.output.reset(TFile::Open(s.candidate.c_str(), "UPDATE"));
  require(s.output && !s.output->IsZombie(), "cannot open copied augmentation candidate");
  s.output->cd();
  s.tree = eventTree(*s.output);
  TLeaf* photons = s.input->GetLeaf("pho_et");
  require(photons && photons->GetLeafCount(), "photon collection counter missing");
  s.photonCount = photons->GetLeafCount();
  for (auto& added : s.added)
  {
    added.buffer.assign(1, 0.);
    s.tree->Branch(("n_" + added.name).c_str(), &added.count, ("n_" + added.name + "/I").c_str());
    s.tree->Branch(added.name.c_str(), added.buffer.data(), (added.name + "[n_" + added.name + "]/D").c_str());
    for (const auto& item : s.request.at("added"))
      if (item.at("name") == added.name && item.contains("title"))
        s.tree->GetBranch(added.name.c_str())->SetTitle(item.at("title").get<std::string>().c_str());
  }
  m_state = std::move(state);
}

BranchAddition::~BranchAddition() = default;

std::string BranchAddition::status() const
{
  require(m_state != nullptr, "branch addition is not open");
  return json{{"entries", m_state->input->GetEntries()}, {"next", m_state->next}}.dump();
}

void BranchAddition::append(const std::string& branch, long long events, const double* values, long long valueCount,
                            const long long* counts)
{
  require(m_state != nullptr, "branch addition is not open");
  State& s = *m_state;
  require(!s.finished, "branch addition already finished");
  require(events >= 0 && s.next + events <= s.input->GetEntries(), "update calculation changed event count");
  auto it = std::find_if(s.added.begin(), s.added.end(), [&](const auto& a) { return a.name == branch; });
  require(it != s.added.end(), "undeclared added branch " + branch);
  require(it->counts.empty(), "values for " + branch + " already supplied for this batch");
  long long total = 0;
  for (long long i = 0; i < events; ++i)
  {
    require(counts[i] >= 0, "negative count for " + branch);
    total += counts[i];
  }
  require(total == valueCount, "value count differs from the counts of " + branch);
  it->counts.assign(counts, counts + events);
  it->values.assign(values, values + total);
}

void BranchAddition::commit(long long events)
{
  require(m_state != nullptr, "branch addition is not open");
  State& s = *m_state;
  require(!s.finished && events >= 0 && s.next + events <= s.input->GetEntries(),
          "update calculation changed event count");
  for (auto& added : s.added)
    require(static_cast<long long>(added.counts.size()) == events,
            "every added branch needs values for the batch: " + added.name);
  for (long long i = 0; i < events; ++i)
  {
    require(s.photonCount->GetBranch()->GetEntry(s.next + i) >= 0, "cannot read photon count");
    const auto photonCount = static_cast<long long>(s.photonCount->GetValue());
    for (auto& added : s.added)
    {
      const long long count = added.counts[static_cast<std::size_t>(i)];
      require(count == photonCount, "score collection is not photon-aligned: " + added.name);
      require(count <= std::numeric_limits<Int_t>::max(), "collection exceeds ROOT counter range: " + added.name);
      if (static_cast<std::size_t>(count) > added.buffer.size())
      {
        added.buffer.resize(static_cast<std::size_t>(count));
        s.tree->SetBranchAddress(added.name.c_str(), added.buffer.data());
      }
      added.count = static_cast<Int_t>(count);
      std::copy(added.values.begin() + static_cast<std::ptrdiff_t>(added.valueOffset),
                added.values.begin() + static_cast<std::ptrdiff_t>(added.valueOffset + count), added.buffer.begin());
      added.valueOffset += static_cast<std::size_t>(count);
      s.output->cd();
      require(s.tree->GetBranch(("n_" + added.name).c_str())->Fill() >= 0 &&
              s.tree->GetBranch(added.name.c_str())->Fill() >= 0, "added branch fill failed");
    }
  }
  for (auto& added : s.added)
  {
    require(added.valueOffset == added.values.size(), "value count differs from the counts of " + added.name);
    added.counts.clear();
    added.values.clear();
    added.valueOffset = 0;
  }
  s.next += events;
}

std::string BranchAddition::finish()
{
  require(m_state != nullptr, "branch addition is not open");
  State& s = *m_state;
  require(!s.finished, "branch addition already finished");
  require(s.next == s.input->GetEntries(), "update omitted input events");
  s.finished = true;
  json metadata = s.request.at("metadata");
  json counters = json::object();
  for (TObject* item : *s.tree->GetListOfBranches())
  {
    auto* branch = static_cast<TBranch*>(item);
    if (branch->GetListOfLeaves()->GetEntries() != 1) continue;
    if (TLeaf* count = static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0))->GetLeafCount())
      counters[branch->GetName()] = count->GetName();
  }
  metadata["schema"]["counters"] = counters;
  {
    NoHistogramDirectory guard;
    s.output->cd();
    require(s.tree->Write("", TObject::kOverwrite) > 0, "TTree::Write failed");
    const json originalMetadata = readMetadata(*s.source);
    for (const auto& [key, value] : metadata.items())
      if (!originalMetadata.contains(key) || originalMetadata.at(key) != value)
        writeText(*s.output, "meta/" + key, dumps(value));
    require(!s.output->TestBit(TFile::kWriteError), "ROOT file write error");
    s.output->Close();
  }
  s.source->Close();
  json added = json::array();
  for (const auto& a : s.added) added.push_back(a.name);
  const std::string proof = provePreserved(
      json{{"original", s.original}, {"candidate", s.candidate}, {"added", added}, {"metadata", metadata}}.dump());
  require(fileSha256(s.original) == s.request.at("input_sha256").get<std::string>(),
          "update input changed during execution");
  return proof;
}

std::string provePreserved(const std::string& requestJson)
{
  const json request = parse(requestJson, "preservation request");
  std::set<std::string> added;
  for (const auto& name : request.at("added")) added.insert(name.get<std::string>());
  auto before = openRead(request.at("original").get<std::string>());
  auto after = openRead(request.at("candidate").get<std::string>());
  require(readMetadata(*after) == request.at("metadata"), "updated metadata differs from the declared change");
  TTree* old = eventTree(*before);
  TTree* now = eventTree(*after);
  require(old->GetEntries() == now->GetEntries(), "update changed event count");
  Columns a(old), b(now);
  for (const auto& buffer : a.all())
  {
    const Buffer* saved = b.find(buffer->name);
    require(saved != nullptr, "update removed an existing branch: " + buffer->name);
    require(branchSignature(*buffer) == branchSignature(*saved), "update changed branch type: " + buffer->name);
  }
  for (const auto& buffer : b.all())
    if (!a.find(buffer->name))
      require(added.count(buffer->name) || (buffer->name.rfind("n_", 0) == 0 && added.count(buffer->name.substr(2))),
              "update added an undeclared branch: " + buffer->name);
  const long long entries = old->GetEntries();
  for (long long entry = 0; entry < entries; ++entry)
  {
    a.read(entry);
    b.read(entry);
    for (const auto& buffer : a.all())
      require(sameValues(*buffer, *b.find(buffer->name)), "update changed protected branch: " + buffer->name);
  }
  // Carried objects other than the direct metadata are byte preserved.
  for (const auto& [name, key] : carriedObjects(*before, false))
  {
    TKey* copy = latestKeyAt(*after, name);
    require(copy != nullptr, "changed carried object: " + name);
    auto x = readObject(key);
    auto y = readObject(copy);
    require(streamed(*x) == streamed(*y), "changed carried object: " + name);
  }
  return json{{"status", "PASS"},
              {"entries", old->GetEntries()},
              {"protected_branches", static_cast<long long>(a.all().size())}}
      .dump();
}

// ===========================================================================
// Group publication
// ===========================================================================
std::string resolveGroup(const std::string& group) { return resolveGroupJson(group).dump(); }

std::string publishGroup(const std::string& requestJson)
{
  const json request = parse(requestJson, "publication request");
  const int stage = request.at("stage").get<int>();
  require(stage == 2 || stage == 3, "invalid publication stage");
  const json& files = request.at("files");
  std::set<std::string> expected, logical, names;
  for (const auto& id : request.at("expected_sources")) expected.insert(id.get<std::string>());
  for (const auto& file : files)
  {
    logical.insert(file.at("logical_id").get<std::string>());
    names.insert(file.at("relative_path").get<std::string>());
  }
  require(files.size() == expected.size() && logical == expected, "group logical-source roster differs");
  require(names.size() == files.size(), "colliding logical file names");
  const fs::path group = fs::absolute(request.at("group").get<std::string>()).lexically_normal();
  fs::create_directories(group);
  struct stat groupInfo{};
  require(::stat(group.c_str(), &groupInfo) == 0, "cannot stat group " + group.string());

  std::vector<json> sorted(files.begin(), files.end());
  std::sort(sorted.begin(), sorted.end(),
            [](const json& a, const json& b) { return a.at("logical_id") < b.at("logical_id"); });
  std::map<std::string, std::array<long long, 4>> identities;
  json records = json::array();
  for (const auto& file : sorted)
  {
    safeRelative(file.at("relative_path").get<std::string>());
    const std::string path = file.at("path").get<std::string>();
    const auto before = fileIdentity(path);
    require(before[0] == static_cast<long long>(groupInfo.st_dev),
            "publication candidate and version root must share a filesystem");
    require(fileSha256(path) == file.at("sha256").get<std::string>(), "group input hash differs");
    require(fileIdentity(path) == before, "group input changed during hashing");
    identities[file.at("logical_id").get<std::string>()] = before;
    require(fileSha256(file.at("receipt_path").get<std::string>()) == file.at("receipt_sha256").get<std::string>(),
            "group validation receipt changed");
    const json receipt = readJsonFile(file.at("receipt_path").get<std::string>());
    require(receipt.value("status", "") == "PASS" &&
                receipt.value("logical_source_id", json()) == file.at("logical_id") &&
                receipt.value("output_sha256", json()) == file.at("sha256"),
            "source validation receipt does not bind file");
    require(receipt.value("proof", json::object()).value("status", "") == "PASS" ||
                receipt.value("proof", json::object()).value("ok", false) == true,
            "source validation receipt has no passing proof");
    records.push_back({{"logical_id", file.at("logical_id")},
                       {"relative_path", file.at("relative_path")},
                       {"sha256", file.at("sha256")},
                       {"receipt_sha256", file.at("receipt_sha256")}});
  }
  const fs::path checks = group / ".validation" / newToken();
  const json population = populationCheck(json(sorted), checks);
  std::error_code ignored;
  fs::remove(checks, ignored);
  json body = {{"schema", "RJDatasetGroupV1"},
               {"group_id", request.at("group_id")},
               {"stage", stage},
               {"recipe", request.at("recipe")},
               {"base_version", request.value("base_version", json())},
               {"files", records},
               {"population", population},
               {"exclusions", request.value("exclusions", json::array())}};
  const std::string version = eventtree::textSha256(dumps(body));
  json manifest = body;
  manifest["version"] = version;
  const fs::path folder = group / "versions" / version;
  fs::create_directories(folder);
  for (const auto& file : sorted)
  {
    const fs::path target = folder / safeRelative(file.at("relative_path").get<std::string>());
    fs::create_directories(target.parent_path());
    const auto& identity = identities.at(file.at("logical_id").get<std::string>());
    if (::link(file.at("path").get<std::string>().c_str(), target.c_str()) != 0)
      require(errno == EEXIST && fileIdentity(target.string()) == identity, "version file collision");
    require(fileIdentity(target.string()) == identity, "version file changed while publishing");
  }
  const std::string raw = dumps(manifest) + "\n";
  const fs::path manifestPath = folder / "MANIFEST.json";
  const fs::path candidate = folder / (".MANIFEST." + newToken());
  {
    const int fd = ::open(candidate.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
    require(fd >= 0, "cannot create manifest candidate");
    const bool written = ::write(fd, raw.data(), raw.size()) == static_cast<ssize_t>(raw.size()) && ::fsync(fd) == 0;
    ::close(fd);
    require(written, "cannot write manifest candidate");
  }
  if (::link(candidate.c_str(), manifestPath.c_str()) != 0)
  {
    require(errno == EEXIST, "cannot publish manifest");
    std::ifstream in(manifestPath, std::ios::binary);
    std::ostringstream existing;
    existing << in.rdbuf();
    require(existing.str() == raw, "manifest collision");
  }
  // One advisory lock per group: independent groups publish concurrently.
  const fs::path lockPath = group / ".publication.lock";
  const int lock = ::open(lockPath.c_str(), O_RDWR | O_CREAT, 0644);
  require(lock >= 0, "cannot open publication lock");
  struct Unlock
  {
    int fd;
    ~Unlock()
    {
      ::flock(fd, LOCK_UN);
      ::close(fd);
    }
  } unlock{lock};
  require(::flock(lock, LOCK_EX) == 0, "cannot lock group publication");
  json previous;
  if (fs::is_symlink(group / "CURRENT")) previous = resolveGroupJson(group.string()).at("manifest");
  const json actual = previous.is_null() ? json() : previous.at("version");
  if (actual == version) return manifest.dump();
  require(actual == request.value("base_version", json()), "publication race: base version changed");
  if (stage == 3)
  {
    require(!previous.is_null(), "augmentation requires a published base");
    json a = json::array(), b = json::array();
    for (const auto& r : previous.at("files")) a.push_back({r.at("logical_id"), r.at("relative_path")});
    for (const auto& r : records) b.push_back({r.at("logical_id"), r.at("relative_path")});
    require(a == b, "augmentation changes logical file set");
  }
  // Recheck the linked version under the lock; candidates are create-once.
  for (const auto& file : sorted)
  {
    const fs::path target = folder / safeRelative(file.at("relative_path").get<std::string>());
    require(fileIdentity(target.string()) == identities.at(file.at("logical_id").get<std::string>()),
            "candidate changed before pointer publication");
    fsyncPath(target.string());
  }
  const fs::path temporary = group / (".CURRENT." + newToken());
  fs::create_directory_symlink(fs::path("versions") / version, temporary);
  fs::rename(temporary, group / "CURRENT");
  fsyncPath(group.string());
  return manifest.dump();
}
}  // namespace eventtree
