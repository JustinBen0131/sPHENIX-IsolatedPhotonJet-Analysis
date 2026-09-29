// eventTree writer core: fills one TTree from flat column buffers.
//
// Each column is a flat array of values plus, for collection branches, entry offsets (entries + 1).
// Scalars and fixed-length arrays are indexed by entry. The Python side owns the buffers and keeps
// them alive until Writer::fill() returns; nothing here copies or reorders values.
#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObjString.h>
#include <TDirectory.h>

#include <algorithm>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <variant>
#include <type_traits>
#include <cmath>

namespace eventtree {

struct Column {
  virtual ~Column() = default;
  virtual void load(long long entry) = 0;
};

template <typename T>
struct Scalar : Column {
  const T* data;
  T value{};
  Scalar(TTree* tree, const std::string& name, const std::string& leaf, const T* d) : data(d) {
    if (auto* b = tree->GetBranch(name.c_str())) b->SetAddress(&value);
    else tree->Branch(name.c_str(), &value, (name + "/" + leaf).c_str());
  }
  void load(long long entry) override { value = data[entry]; }
};

template <typename T>
struct Fixed : Column {
  const T* data;
  int size;
  std::vector<T> value;
  Fixed(TTree* tree, const std::string& name, const std::string& leaf, const T* d, int n)
      : data(d), size(n), value(n) {
    if (auto* b = tree->GetBranch(name.c_str())) b->SetAddress(value.data());
    else tree->Branch(name.c_str(), value.data(), (name + "[" + std::to_string(n) + "]/" + leaf).c_str());
  }
  void load(long long entry) override {
    std::memcpy(value.data(), data + entry * size, size * sizeof(T));
  }
};

// A collection branch is a C array whose length is an Int_t counter branch shared by its group. An
// empty collection then costs almost nothing per entry, unlike a std::vector's per-entry header.
struct Counter : Column {
  const long long* offsets;
  Int_t value = 0;
  Counter(TTree* tree, const std::string& name, const long long* o) : offsets(o) {
    if (auto* b = tree->GetBranch(name.c_str())) b->SetAddress(&value);
    else tree->Branch(name.c_str(), &value, (name + "/I").c_str());
  }
  void load(long long entry) override { value = static_cast<Int_t>(offsets[entry + 1] - offsets[entry]); }
};

template <typename T>
struct Counted : Column {
  const T* data;
  const long long* offsets;
  std::vector<T> buffer;
  Counted(TTree* tree, const std::string& name, const std::string& leaf, const T* d, const long long* o,
          const std::string& counter, long long max_count)
      : data(d), offsets(o), buffer(max_count > 0 ? max_count : 1) {
    if (auto* b = tree->GetBranch(name.c_str())) b->SetAddress(buffer.data());
    else tree->Branch(name.c_str(), buffer.data(), (name + "[" + counter + "]/" + leaf).c_str());
  }
  void load(long long entry) override {
    std::copy(data + offsets[entry], data + offsets[entry + 1], buffer.begin());
  }
};

struct StringVector : Column {
  std::vector<std::string> data;
  std::vector<long long> offsets;
  std::vector<std::string> value;
  std::vector<std::string>* address = &value;
  StringVector(TTree* tree, const std::string& name, std::vector<std::string> d, std::vector<long long> o)
      : data(std::move(d)), offsets(std::move(o)) {
    if (auto* b = tree->GetBranch(name.c_str())) b->SetAddress(&address);
    else tree->Branch(name.c_str(), &address);
  }
  void load(long long entry) override {
    value.assign(data.begin() + offsets[entry], data.begin() + offsets[entry + 1]);
  }
};

class Writer {
 public:
  Writer(const std::string& path, const std::string& tree_name, const std::string& title, int compression,
         long long entries)
      : entries_(entries) {
    file_.reset(TFile::Open(path.c_str(), "CREATE"));
    if (!file_ || file_->IsZombie()) throw std::runtime_error("cannot create " + path);
    file_->SetCompressionSettings(compression);
    tree_ = new TTree(tree_name.c_str(), title.c_str());
    tree_->SetAutoFlush(-16 * 1024 * 1024);
  }

  // All previous fill calls finished before Python releases their buffers.
  // Existing branches are rebound, so only one bounded batch is held at once.
  void beginBatch(long long entries) {
    tree_->ResetBranchAddresses();
    columns_.clear();
    entries_ = entries;
  }

  template <typename T>
  void scalar(const std::string& name, const std::string& leaf, const T* data, const std::string& title) {
    columns_.emplace_back(new Scalar<T>(tree_, name, leaf, data));
    retitle(name, title);
  }
  template <typename T>
  void fixed(const std::string& name, const std::string& leaf, const T* data, int size, const std::string& title) {
    columns_.emplace_back(new Fixed<T>(tree_, name, leaf, data, size));
    retitle(name, title);
  }
  void counter(const std::string& name, const long long* offsets) {
    columns_.emplace_back(new Counter(tree_, name, offsets));
  }
  template <typename T>
  void counted(const std::string& name, const std::string& leaf, const T* data, const long long* offsets,
               const std::string& counter, long long max_count, const std::string& title) {
    columns_.emplace_back(new Counted<T>(tree_, name, leaf, data, offsets, counter, max_count));
    retitle(name, title);
  }
  void strings(const std::string& name, std::vector<std::string> data, std::vector<long long> offsets,
               const std::string& title) {
    columns_.emplace_back(new StringVector(tree_, name, std::move(data), std::move(offsets)));
    retitle(name, title);
  }

  long long fill() {
    for (long long entry = 0; entry < entries_; ++entry) {
      for (auto& column : columns_) column->load(entry);
      if (tree_->Fill() < 0) throw std::runtime_error("TTree::Fill failed");
    }
    return entries_;
  }

  // Metadata is written into ordinary ROOT objects, never extra event rows.
  void text(const std::string& path, const std::string& value) {
    TObjString text(value.c_str());
    object(path,text);
  }
  void object(const std::string& path, TObject& value) {
    if(path.empty() || path==tree_->GetName() || path.find("..")!=std::string::npos)
      throw std::runtime_error("invalid metadata object path: "+path);
    const auto separator = path.find_last_of('/');
    const std::string directory = separator == std::string::npos ? "" : path.substr(0, separator);
    const std::string name = separator == std::string::npos ? path : path.substr(separator + 1);
    TDirectory* target = file_.get();
    if (!directory.empty()) {
      target = file_->GetDirectory(directory.c_str());
      if (!target) {
        file_->mkdir(directory.c_str(),"",true);
        target = file_->GetDirectory(directory.c_str());
      }
    }
    if (!target) throw std::runtime_error("cannot create metadata directory: " + directory);
    target->cd();
    if (value.Write(name.c_str(), TObject::kOverwrite) <= 0)
      throw std::runtime_error("cannot write metadata: " + path);
    file_->cd();
  }

  void close() {
    file_->cd();
    // kOverwrite replaces the cycle an autosave wrote during fill(): the file holds exactly one eventTree key
    if (tree_->Write("", TObject::kOverwrite) <= 0) throw std::runtime_error("TTree::Write failed");
    if (file_->TestBit(TFile::kWriteError)) throw std::runtime_error("ROOT file write error");
    file_->Close();
    if (file_->TestBit(TFile::kWriteError)) throw std::runtime_error("ROOT file close error");
  }

 private:
  void retitle(const std::string& name, const std::string& title) {
    if (!title.empty()) tree_->GetBranch(name.c_str())->SetTitle(title.c_str());
  }
  std::unique_ptr<TFile> file_;
  TTree* tree_ = nullptr;
  std::vector<std::unique_ptr<Column>> columns_;
  long long entries_;
};

// A single event's typed columns. Production feeds its in-memory records into
// this same writer core used by Python consolidation and augmentation. Schema
// declarations determine storage width; integer identities never pass through
// floating point. Empty collections keep their declared type.
struct FieldSpec {
  std::string name, type, title, counter;
  int length = 0;  // 0 scalar, -1 counted vector, positive fixed array
};

using FieldValues = std::variant<std::vector<std::int8_t>, std::vector<std::uint8_t>,
    std::vector<std::int16_t>, std::vector<std::uint16_t>, std::vector<std::int32_t>,
    std::vector<std::uint32_t>, std::vector<Long64_t>, std::vector<ULong64_t>,
    std::vector<float>, std::vector<double>, std::vector<std::string>>;

struct Field {
  FieldSpec spec;
  FieldValues values;
  std::array<long long, 2> offsets{{0, 0}};
  bool assigned = false;
  std::string leaf;

  explicit Field(FieldSpec s) : spec(std::move(s)) {
    if (spec.type == "int8") { values = std::vector<std::int8_t>(); leaf = "B"; }
    else if (spec.type == "uint8") { values = std::vector<std::uint8_t>(); leaf = "b"; }
    else if (spec.type == "int16") { values = std::vector<std::int16_t>(); leaf = "S"; }
    else if (spec.type == "uint16") { values = std::vector<std::uint16_t>(); leaf = "s"; }
    else if (spec.type == "int32") { values = std::vector<std::int32_t>(); leaf = "I"; }
    else if (spec.type == "uint32") { values = std::vector<std::uint32_t>(); leaf = "i"; }
    else if (spec.type == "int64") { values = std::vector<Long64_t>(); leaf = "L"; }
    else if (spec.type == "uint64") { values = std::vector<ULong64_t>(); leaf = "l"; }
    else if (spec.type == "float") { values = std::vector<float>(); leaf = "F"; }
    else if (spec.type == "double") { values = std::vector<double>(); leaf = "D"; }
    else if (spec.type == "string" && spec.length == -1) values = std::vector<std::string>();
    else throw std::runtime_error("unsupported eventTree field type: " + spec.type);
  }

  template<class T> void append(T value) {
    std::visit([&](auto& data) {
      using Value = typename std::decay_t<decltype(data)>::value_type;
      if constexpr (std::is_arithmetic_v<T> && std::is_arithmetic_v<Value>) {
        if constexpr (std::is_integral_v<Value>) {
          bool invalid = false;
          if constexpr (std::is_integral_v<T>) {
            if constexpr (std::is_signed_v<T>) {
              if (value < 0) {
                if constexpr (!std::is_signed_v<Value>) invalid = true;
                else invalid = static_cast<std::intmax_t>(value) < std::numeric_limits<Value>::lowest();
              } else invalid = static_cast<std::uintmax_t>(value) > static_cast<std::uintmax_t>(std::numeric_limits<Value>::max());
            } else invalid = static_cast<std::uintmax_t>(value) > static_cast<std::uintmax_t>(std::numeric_limits<Value>::max());
          } else if constexpr (std::is_floating_point_v<T>) {
            const long double wide = static_cast<long double>(value);
            const long double upper = std::ldexp(1.0L, std::numeric_limits<Value>::digits);
            const long double lower = std::is_signed_v<Value> ? -upper : 0.0L;
            invalid = !std::isfinite(wide) || std::trunc(wide) != wide || wide < lower || wide >= upper;
          }
          if (invalid) throw std::runtime_error("integer overflow in " + spec.name);
        }
        data.push_back(static_cast<Value>(value));
      } else if constexpr (std::is_same_v<T, std::string> && std::is_same_v<Value, std::string>) data.push_back(value);
      else throw std::runtime_error("eventTree field type mismatch: " + spec.name);
    }, values);
    assigned = true;
  }

  void reset() { std::visit([](auto& data) { data.clear(); }, values); assigned = false; }

  void bind(Writer& writer) {
    std::visit([&](auto& data) {
      using Value = typename std::decay_t<decltype(data)>::value_type;
      if (spec.length >= 0 && data.size() != static_cast<std::size_t>(std::max(1, spec.length)))
        throw std::runtime_error("scalar/fixed field length differs: " + spec.name);
      offsets[1] = static_cast<long long>(data.size());
      if constexpr (std::is_same_v<Value, std::string>) {
        writer.strings(spec.name, data, {0, offsets[1]}, spec.title);
      } else {
        static const Value empty{};
        const Value* pointer = data.empty() ? &empty : data.data();
        if (spec.length == -1) {
          if (data.size() > static_cast<std::size_t>(std::numeric_limits<Int_t>::max()))
            throw std::runtime_error("collection exceeds ROOT counter range: " + spec.name);
          writer.counter(spec.counter, offsets.data());
          writer.counted(spec.name, leaf, pointer, offsets.data(), spec.counter, offsets[1], spec.title);
        } else if (spec.length == 0) writer.scalar(spec.name, leaf, pointer, spec.title);
        else writer.fixed(spec.name, leaf, pointer, spec.length, spec.title);
      }
    }, values);
  }
};

class Entry {
 public:
  explicit Entry(const std::vector<FieldSpec>& specs) {
    for (const auto& spec : specs) {
      if (!fields_.emplace(spec.name, Field(spec)).second)
        throw std::runtime_error("duplicate eventTree field: " + spec.name);
      order_.push_back(spec.name);
    }
  }
  void reset() { for (auto& item : fields_) item.second.reset(); }
  template<class T> void put(const std::string& name, T value) { fields_.at(name).append(value); }
  template<class Values> void putAll(const std::string& name, const Values& values) {
    auto& field = fields_.at(name); field.assigned = true;
    for (const auto& value : values) field.append(value);
  }
  void empty(const std::string& name) {
    auto& field = fields_.at(name);
    if (field.spec.length != -1) throw std::runtime_error("only vector fields may be declared empty: " + name);
    field.assigned = true;
  }
  std::vector<std::string> missing() const {
    std::vector<std::string> names;
    for (const auto& item : fields_) if (!item.second.assigned) names.push_back(item.first);
    return names;
  }
  // Every field must have an explicit filling rule, including empty vectors.
  // This guard prevents default zero/NaN buffers from masquerading as capture.
  void write(Writer& writer, bool bookOnly = false) {
    const auto gaps = missing();
    if (!gaps.empty()) throw std::runtime_error("unfilled eventTree field: " + gaps.front());
    writer.beginBatch(bookOnly ? 0 : 1);
    for (const auto& name : order_) fields_.at(name).bind(writer);
    writer.fill();
  }
 private:
  std::map<std::string, Field> fields_;
  std::vector<std::string> order_;
};

}  // namespace eventtree
