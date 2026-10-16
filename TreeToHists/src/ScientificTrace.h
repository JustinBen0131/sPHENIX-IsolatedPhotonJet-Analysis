/* Scientific values and effective definitions shared by the native producers.
   This header performs no selection, weighting, fitting or I/O. Its arguments
   are the values used by the producer. File/run/source digests are provenance,
   never substitutes for a scientific definition. */
#ifndef PHOTONJET_SCIENTIFIC_TRACE_H
#define PHOTONJET_SCIENTIFIC_TRACE_H
#include <functional>
#include <nlohmann/json.hpp>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace photonjet_science {
using Json = nlohmann::json;
inline const std::vector<std::string>& categories() {
  static const std::vector<std::string> names = {"selection", "weights", "models", "normalization", "projection"};
  return names;
}
inline std::string real(double x) {
  static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559, "IEEE754 binary64 is required");
  if (!std::isfinite(x)) throw std::runtime_error("nonfinite scientific value requires an explicit measurement state");
  if (x==0) x=0;
  std::uint64_t bits; std::memcpy(&bits,&x,8);
  std::ostringstream s; s << "f64:" << std::hex << std::setw(16) << std::setfill('0') << bits; return s.str();
}
inline Json portable(const Json& value) {
  if (value.is_number_float()) return real(value.get<double>());
  Json result=value;
  if (result.is_structured()) for (auto& child:result) child=portable(child);
  return result;
}
inline Json numbers(const std::vector<double>& values) {
  Json a=Json::array(); for(double value:values)a.push_back(real(value)); return a;
}
inline void require(bool condition,const std::string& reason) {
  if(!condition)throw std::runtime_error("ScientificTrace: "+reason);
}
inline bool containsUnknown(const Json& value) {
  if(value.is_object() && (value.value("state",std::string())=="unknown" || value.value("state",std::string())=="missing"))return true;
  if(value.is_structured())for(const auto& v:value)if(containsUnknown(v))return true;
  return false;
}
inline Json unknown(const std::string& reason) {
  require(!reason.empty(),"unknown needs a reason"); return {{"state","unknown"},{"reason",reason}};
}
inline Json known(const Json& value) { return {{"state","measured"},{"value",portable(value)}}; }
inline Json absent(const std::string& reason) {
  require(!reason.empty(),"not applicable needs an explicit reason");
  return {{"state","not_applicable"},{"reason",reason}};
}

/* Terms carry exactly the configuration or model component consumed by an
   observable. A term applies only where its axis windows intersect the bins
   consumed by that observable. Windows are half-open physical coordinates;
   no histogram/file name is parsed to invent a cut. */
inline bool intersects(const Json& termWindows,const Json& selected) {
  for(const auto& a:termWindows)for(const auto& b:selected)if(a.at("axis")==b.at("axis")) {
    const double low=std::max(a.at("low").get<double>(),b.at("low").get<double>());
    const double high=std::min(a.at("high").get<double>(),b.at("high").get<double>());
    if(!(low<high))return false;
  }
  return true;
}
inline Json window(int axis,double low,double high) { return Json::array({{{"axis",axis},{"low",low},{"high",high}}}); }
inline Json term(const std::string& category,const std::string& id,const Json& value,
                 const Json& windows=Json::array(),const std::string& reference="",const Json& origin=Json()) {
  require(!id.empty(),"empty term identity");
  require(std::find(categories().begin(),categories().end(),category)!=categories().end(),"unknown term category");
  require(windows.is_array(),"term windows must be an array");
  for(const auto& w:windows)require(w.at("axis").is_number_integer() && w.at("axis").get<int>()>=0 &&
    w.at("axis").get<int>()<3 && w.at("low").get<double>()<w.at("high").get<double>(),"invalid term window");
  Json out={{"category",category},{"id",id},{"value",portable(value)},{"windows",windows}};
  if(!reference.empty())out["reference"]=reference;
  if(!origin.is_null()) {require(!reference.empty(),"origin requires a stable definition reference");out["origin"]=origin;}
  return out;
}
inline Json definition(const std::string& observable,const Json& identity,const Json& terms,
                       const Json& declarations) {
  require(!observable.empty() && identity.is_object() && !identity.empty(),"observable and scientific channel required");
  require(terms.is_array() && declarations.is_object(),"definition shape");
  std::set<std::string> ids;
  for(const auto& item:terms)require(ids.insert(item.at("category").get<std::string>()+"/"+item.at("id").get<std::string>()).second,
                                  "duplicate definition term");
  for(const auto& c:categories())require(declarations.contains(c),"missing category declaration "+c);
  return {{"schema","EffectiveObservableDefinitionV1"},{"observable",observable},{"channel",portable(identity)},
          {"terms",terms},{"categories",declarations}};
}

/* The native trace is attached to a value, not kept in global last-read state.
   A projection filters the definition's terms before they enter physics
   identity. Source file paths and hashes may accompany provenance elsewhere. */
struct Trace {
  Json roots=Json::object(), operations=Json::array(), missing=Json::array(), origins=Json::object();
  static Trace unavailable(const std::string& why) { Trace t;t.missing.push_back(why);return t; }
  static Trace source(const Json& d,const Json& windows=Json::array()) {
    if(d.is_null())return unavailable("producer did not emit an effective observable definition");
    if(d.value("schema","")=="EffectiveObservableDefinitionsV1") {
      const auto& definitions=d.at("definitions");
      if(!definitions.is_array()||definitions.empty())return unavailable("record population has no observed definitions");
      std::vector<Trace> inputs;for(const auto& item:definitions)inputs.push_back(source(item,windows));
      return combine("union_record_populations",inputs);
    }
    require(d.value("schema","")=="EffectiveObservableDefinitionV1","wrong producer definition schema");
    Trace t;
    Json selected=Json::object(),references=Json::object();
    for(const auto& entry:d.at("terms"))if(intersects(entry.at("windows"),windows)) {
      const std::string key=entry.at("category").get<std::string>()+"/"+entry.at("id").get<std::string>();
      selected[key]=entry.at("value");
      if(entry.contains("reference")) {
        const std::string reference=entry.at("reference");references[key]=reference;
        if(entry.contains("origin"))t.origins[reference]=entry.at("origin");
      }
    }
    const Json key={{"observable",d.at("observable")},{"channel",d.at("channel")}};
    t.roots[key.dump()]={{"identity",key},{"terms",selected},{"categories",d.at("categories")},
                         {"consumed_windows",portable(windows)},{"references",references}};
    return t;
  }
  bool complete(const std::string& category) const {
    if(!missing.empty() || roots.empty())return false;
    for(const auto& root:roots) {
      const auto row=root.at("categories").value(category,unknown("missing declaration"));
      const auto state=row.value("state","");
      if(state!="measured" && state!="not_applicable")return false;
      for(const auto& [key,value]:root.at("terms").items())if(key.rfind(category+"/",0)==0 && containsUnknown(value))return false;
    }
    return true;
  }
  Trace then(const std::string& operation,const Json& parameters=Json::object()) const {
    require(!operation.empty(),"operation identity is required");
    Trace out=*this;out.operations.push_back({{"id",operation},{"parameters",portable(parameters)}});return out;
  }
  /* A consumed canonical calculation setting is an input edge, not merely
     descriptive operation text. The stable reference excludes its value.
     Origin describes where the native value was read; qualification binds it
     to the exact source/configuration span. Omitted origin never authorizes a
     canonical edit, although the measured value still affects comparison. */
  Trace withParameter(const std::string& reference,const std::string& category,
                      const Json& value,const Json& origin=Json()) const {
    require(!reference.empty(),"calculation parameter requires a stable reference");
    Json declarations=Json::object();
    for(const auto& c:categories())declarations[c]=c==category?known(true):absent("this parameter belongs to "+category);
    const auto parameter=Trace::source(definition("canonical_calculation_parameter",{{"reference",reference}},
      Json::array({term(category,"effective_value",value,Json::array(),reference,origin)}),declarations));
    Trace out=*this;
    for(const auto& [key,row]:parameter.roots.items()) {
      require(!out.roots.contains(key)||out.roots.at(key)==row,"conflicting calculation parameter");out.roots[key]=row;
    }
    for(const auto& [id,site]:parameter.origins.items()) {
      require(!out.origins.contains(id)||out.origins.at(id)==site,"conflicting parameter origin");out.origins[id]=site;
    }
    return out;
  }
  /* Distinct scientific operands can consume different values of the same
     canonical setting. Namespace by the caller's stable semantic role, never
     by the input path/hash or the setting value. Repeated conflicting values
     within one role are still rejected by combine. */
  Trace inRole(const std::string& role) const {
    require(!role.empty(),"trace operand role is empty");
    const auto reference=[&](const std::string& id){return Json({{"role",role},{"reference",id}}).dump();};
    const auto scopedRoots=[&](const Json& inputs) {
      Json result=Json::object();
      for(const auto& row:inputs) {
        Json value=row;auto& identity=value["identity"];
        Json chain=identity.value("roles",Json::array());chain.insert(chain.begin(),role);identity["roles"]=chain;
        for(auto& id:value["references"])id=reference(id.get<std::string>());
        const std::string key=identity.dump()+"/"+value.at("consumed_windows").dump();
        require(!result.contains(key)||result.at(key)==value,"conflicting roots within semantic role");result[key]=value;
      }
      return result;
    };
    std::function<Json(const Json&)> scopedTrace;
    scopedTrace=[&](const Json& input) {
      Json value=input;value["roots"]=scopedRoots(input.at("roots"));
      for(auto& operation:value["operations"])if(operation.contains("operands"))
        for(auto& operand:operation["operands"])operand=scopedTrace(operand);
      return value;
    };
    Trace out=*this;const auto scoped=scopedTrace(scientificJson());
    out.roots=scopedRoots(roots);out.operations=scoped.at("operations");out.origins=Json::object();
    for(const auto& [id,origin]:origins.items())out.origins[reference(id)]=origin;
    return out;
  }
  static Trace combine(const std::string& operation,const std::vector<Trace>& inputs,
                       const Json& parameters=Json::object()) {
    require(!inputs.empty(),"operation has no inputs");
    Trace out;
    Json operands=Json::array();
    for(const auto& input:inputs) {
      operands.push_back(input.scientificJson());
      for(const auto& [id,origin]:input.origins.items()) {
        require(!out.origins.contains(id)||out.origins.at(id)==origin,"conflicting canonical origins");out.origins[id]=origin;
      }
      for(const auto& [key,row]:input.roots.items()) {
        /* Distinct restrictions of one histogram are distinct dependency edges. */
        const std::string edge=key+"/"+row.at("consumed_windows").dump();
        if(out.roots.contains(edge))require(out.roots.at(edge)==row,"contradictory definitions for one dependency");
        else out.roots[edge]=row;
      }
      for(const auto& why:input.missing)out.missing.push_back(why);
    }
    out.operations.push_back({{"id",operation},{"parameters",portable(parameters)},{"operands",operands}});
    return out;
  }
  Json scientificJson() const {
    Json inputs=roots;
    for(auto& root:inputs)root.erase("categories");
    return {{"roots",inputs},{"operations",operations}};
  }
  Json json() const { return {{"roots",roots},{"operations",operations},{"missing",missing},{"origins",origins}}; }
  static Trace fromJson(const Json& j) {
    if(j.is_null())return unavailable("value has no scientific trace");
    Trace t;t.roots=j.at("roots");t.operations=j.at("operations");t.missing=j.at("missing");t.origins=j.value("origins",Json::object());return t;
  }
  Json definitions() const {
    Json selected=Json::object();
    for(const auto& [id,row]:roots.items())selected[id]={{"identity",row.at("identity")},{"terms",row.at("terms")},
                                                       {"consumed_windows",row.at("consumed_windows")}};
    return {{"inputs",selected},{"operations",operations}};
  }
  std::string identity() const {
    Json identities=Json::array();for(const auto& root:roots)identities.push_back(root.at("identity"));
    /* A stable scientific series key excludes actual values, cut parameter
       values and style. A role at the consuming expression disambiguates two
       operations on the same inputs. */
    return identities.dump();
  }
};

inline Json covariance(const std::string& state,const std::vector<double>& data,std::size_t n,
                       const std::string& method) {
  require(!method.empty(),"covariance method/absence reason is required");
  if(state=="measured") {
    require(data.size()==n*n,"covariance shape mismatch");
    for(std::size_t i=0;i<n;++i) {
      require(data[i*n+i]>=0,"negative diagonal variance");
      for(std::size_t j=0;j<n;++j)require(data[i*n+j]==data[j*n+i],"nonsymmetric covariance");
    }
    return {{"state",state},{"dimension",n},{"values",numbers(data)},{"method",method}};
  }
  require(data.empty(),"unmeasured covariance supplied values");
  require(state=="not_computed" || state=="not_applicable" || state=="unknown","covariance availability state");
  return {{"state",state},{"method",method}};
}

/* measured, undefined (e.g. zero denominator), and missing are different.
   Missing is never converted to zero or silently omitted from the roster. */
inline Json bins(const std::vector<double>& edges,const std::vector<double>& values,
                 const std::vector<double>& variance,const std::vector<std::string>& states,
                 const Json& covarianceRecord) {
  require(edges.size()==values.size()+1 && variance.size()==values.size() && states.size()==values.size(),"bin shape mismatch");
  for(std::size_t i=1;i<edges.size();++i)require(edges[i]>edges[i-1],"unordered bin edges");
  Json y=Json::array(),e=Json::array();
  for(std::size_t i=0;i<values.size();++i) {
    require(states[i]=="measured" || states[i]=="undefined" || states[i]=="missing","unknown measurement state");
    if(states[i]=="measured") {
      require(std::isfinite(variance[i]) && variance[i]>=0,"invalid variance");
      y.push_back(real(values[i]));e.push_back(real(variance[i]));
    } else { y.push_back(nullptr);e.push_back(nullptr); }
  }
  return {{"type","binned"},{"edges",numbers(edges)},{"values",y},{"variances",e},{"measurement_state",states},
          {"covariance",covarianceRecord}};
}
inline bool measurementsComplete(const Json& q) {
  if(q.is_object() && q.contains("type")) {
    const std::string type=q.at("type");
    if(type=="binned" && q.at("values").empty())return false;
    if((type=="histogram"||type=="selected_histogram") && q.at("cells").empty())return false;
    if((type=="points"||type=="symmetric_points"||type=="asymmetric_points") && q.at("rows").empty())return false;
  }
  if(q.contains("measurement_state"))for(const auto& state:q.at("measurement_state"))if(state=="missing")return false;
  if(q.contains("state") && (q.at("state")=="missing" || q.at("state")=="unknown"))return false;
  if(q.is_structured())for(const auto& child:q)if(child.is_structured() && !measurementsComplete(child))return false;
  return true;
}
}  // namespace photonjet_science
#endif
