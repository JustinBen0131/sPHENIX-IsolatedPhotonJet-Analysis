/* ROOT adapters for ScientificTrace. Only explicitly bound scientific objects
   enter a complete receipt. Unbound graphs/maps are reported, never guessed
   from their title, colour, ROOT counter or position in the canvas. */
#ifndef PHOTONJET_ROOT_SCIENTIFIC_TRACE_H
#define PHOTONJET_ROOT_SCIENTIFIC_TRACE_H
#include "PlotPhysics.h"
#include "../../TreeToHists/src/ScientificTrace.h"
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TList.h>
#include <TNamed.h>
#include <TPad.h>
#include <TFrame.h>
#include <TObjString.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <set>

namespace photonjet_science {
inline TList* tags(TObject* object) {
  if(auto* h=dynamic_cast<TH1*>(object))return h->GetListOfFunctions();
  if(auto* g=dynamic_cast<TGraph*>(object))return g->GetListOfFunctions();
  if(auto* f=dynamic_cast<TF1*>(object))return f->GetHistogram()->GetListOfFunctions();
  if(auto* t=dynamic_cast<TTree*>(object))return t->GetUserInfo();
  return nullptr;
}
inline Json tag(const TObject* object,const char* name) {
  auto* list=tags(const_cast<TObject*>(object));
  auto* saved=list?dynamic_cast<TNamed*>(list->FindObject(name)):nullptr;
  return saved?Json::parse(saved->GetTitle()):Json();
}
inline void setTag(TObject* object,const char* name,const Json& value) {
  auto* list=tags(object);require(list,"object cannot retain a scientific binding");
  if(auto* old=list->FindObject(name))list->Remove(old);
  list->Add(new TNamed(name,value.dump().c_str()));
}
inline Json axis(const TAxis* a) {
  Json x=Json::array();for(int i=1;i<=a->GetNbins()+1;++i)x.push_back(real(a->GetBinLowEdge(i)));return x;
}
inline Json histogram(const TH1* h,bool flow=false) {
  Json axes=Json::array({axis(h->GetXaxis())});
  if(h->GetDimension()>1)axes.push_back(axis(h->GetYaxis()));
  if(h->GetDimension()>2)axes.push_back(axis(h->GetZaxis()));
  Json rows=Json::array();
  for(int i=0;i<h->GetNcells();++i) {
    if(!flow && (h->IsBinUnderflow(i)||h->IsBinOverflow(i)))continue;
    int x,y,z;h->GetBinXYZ(i,x,y,z);
    const double value=h->GetBinContent(i),error=h->GetBinError(i);
    const bool valid=std::isfinite(value)&&std::isfinite(error)&&error>=0;
    rows.push_back({{"index",{x,y,z}},{"state",valid?"measured":"missing"},
      {"value",valid?Json(real(value)):Json()},{"variance",valid?Json(real(error*error)):Json()}});
  }
  return {{"type","histogram"},{"axes",axes},{"flow",flow?"included":"excluded"},{"cells",rows}};
}
inline Json graph(const TGraph* g) {
  Json rows=Json::array();
  const auto* sym=dynamic_cast<const TGraphErrors*>(g);
  const auto* asym=dynamic_cast<const TGraphAsymmErrors*>(g);
  for(int i=0;i<g->GetN();++i) {
    double x,y;g->GetPoint(i,x,y);
    Json row={{"x",real(x)},{"y",real(y)}};
    if(asym)row["errors"]={real(asym->GetErrorXlow(i)),real(asym->GetErrorXhigh(i)),real(asym->GetErrorYlow(i)),real(asym->GetErrorYhigh(i))};
    else if(sym)row["errors"]={real(sym->GetErrorX(i)),real(sym->GetErrorY(i))};
    rows.push_back(row);
  }
  return {{"type",asym?"asymmetric_points":sym?"symmetric_points":"points"},{"rows",rows}};
}
inline Json nativeValues(const TObject* object) {
  if(const auto* h=dynamic_cast<const TH1*>(object))return histogram(h);
  if(const auto* g=dynamic_cast<const TGraph*>(object))return graph(g);
  if(const auto* f=dynamic_cast<const TF1*>(object)) {
    Json parameters=Json::array();
    for(int i=0;i<f->GetNpar();++i)parameters.push_back({{"slot",i},{"value",real(f->GetParameter(i))},{"error",real(f->GetParError(i))}});
    return {{"type","analytic_function"},{"formula",f->GetExpFormula().Data()},{"range",{real(f->GetXmin()),real(f->GetXmax())}},
      {"parameters",parameters},{"chi2",real(f->GetChisquare())},{"ndf",f->GetNDF()}};
  }
  return Json();
}
inline Json definitionOf(const TObject* object) { return tag(object,"physics_effective_definition"); }
inline Trace traceOf(const TObject* object) {
  const auto derived=tag(object,"physics_trace");
  return derived.is_null()?Trace::source(definitionOf(object)):Trace::fromJson(derived);
}
inline void attachDefinition(TObject* object,const Json& definition) { setTag(object,"physics_effective_definition",definition); }
inline void attachTrace(TObject* object,const Trace& trace) { setTag(object,"physics_trace",trace.json()); }
/* JSON/TParameter producers retain a neighbouring physics_<name> record.
   Its absence stays unknown; the value itself is never an ancestry witness. */
inline Json companion(TDirectory& directory,const std::string& path) {
  const auto slash=path.rfind('/');
  const std::string name=slash==std::string::npos?"physics_"+path:path.substr(0,slash+1)+"physics_"+path.substr(slash+1);
  auto* saved=dynamic_cast<TObjString*>(directory.Get(name.c_str()));
  return saved?Json::parse(saved->GetString().Data()):Json();
}
inline Trace companionTrace(TDirectory& directory,const std::string& path) {
  const auto record=companion(directory,path);
  return record.is_object() && record.contains("trace")?Trace::fromJson(record.at("trace")):
      Trace::unavailable("producer did not retain a trace for "+path);
}
inline Trace projectionOf(const TObject* object,const Json& windows) {
  const auto derived=tag(object,"physics_trace");
  if(!derived.is_null())return Trace::fromJson(derived).then("restrict",{{"windows",windows}});
  const auto d=definitionOf(object);
  if(!d.is_null())return Trace::source(d,windows);
  const auto definitions=tag(object,"physics_effective_definitions");
  if(definitions.is_array() && !definitions.empty()) {
    std::vector<Trace> inputs;for(const auto& definition:definitions)inputs.push_back(Trace::source(definition,windows));
    return Trace::combine("sum_input_populations",inputs);
  }
  const auto constituents=tag(object,"physics_constituents");
  if(!constituents.is_null()) {
    std::vector<Trace> inputs;
    for(const auto& row:constituents) {
      Trace input;
      if(row.contains("definitions") && row.at("definitions").is_array() && !row.at("definitions").empty()) {
        std::vector<Trace> members;for(const auto& definition:row.at("definitions"))members.push_back(Trace::source(definition,windows));
        input=Trace::combine("sum_input_populations",members);
      } else input=Trace::source(row.at("definition"),windows);
      for(const auto& factor:row.at("factors"))input=input.then("sample_normalization",{{"factor",factor}});
      inputs.push_back(input);
    }
    return Trace::combine("sum_independent_sources",inputs);
  }
  return traceOf(object).then("restrict",{{"windows",windows}});
}
inline Trace recordTrace(TDirectory& directory,const std::string& path,const std::string& column="") {
  auto* object=directory.Get(path.c_str());
  if(!column.empty()) {
    const auto traces=tag(object,"physics_column_traces");
    if(traces.is_object() && traces.contains(column))return Trace::fromJson(traces.at(column));
    return Trace::unavailable("producer did not retain a column-specific trace for "+path+"/"+column);
  }
  if(object && tags(object))return projectionOf(object,Json::array());
  return companionTrace(directory,path);
}
inline void sumDefinition(TH1* out,const TH1* input,double factor,bool first) {
  Json rows=first?Json::array():tag(out,"physics_constituents");
  if(rows.is_null())rows=Json::array();
  auto incoming=tag(input,"physics_constituents");
  if(incoming.is_null())rows.push_back({{"definition",definitionOf(input)},{"definitions",tag(input,"physics_effective_definitions")},{"factors",Json::array({real(factor)})}});
  else for(auto row:incoming) { row["factors"].push_back(real(factor));rows.push_back(row); }
  if(auto* old=tags(out)->FindObject("physics_effective_definition"))tags(out)->Remove(old);
  if(auto* old=tags(out)->FindObject("physics_effective_definitions"))tags(out)->Remove(old);
  if(auto* old=tags(out)->FindObject("physics_trace"))tags(out)->Remove(old);
  setTag(out,"physics_constituents",rows);
}

/* A renderer consumes a specific physical interval. Values outside that
   interval, marker offsets and extra points closing a filled shape are absent
   from its scientific series. The calculation's full covariance remains in
   the bound upstream tag when available. */
inline Json selectedHistogram(const TH1* h,double low,double high) {
  require(h && h->GetDimension()==1,"selected histogram must be one dimensional");
  Json cells=Json::array();
  for(int i=1;i<=h->GetNbinsX();++i) {
    const double a=h->GetXaxis()->GetBinLowEdge(i),b=h->GetXaxis()->GetBinUpEdge(i);
    if(a<low-1e-9 || b>high+1e-9)continue;
    const double y=h->GetBinContent(i),e=h->GetBinError(i);
    const bool valid=std::isfinite(y)&&std::isfinite(e)&&e>=0;
    cells.push_back({{"low",real(a)},{"high",real(b)},{"state",valid?"measured":"missing"},
      {"value",valid?Json(real(y)):Json()},{"variance",valid?Json(real(e*e)):Json()}});
  }
  return {{"type","selected_histogram"},{"cells",cells}};
}
inline Json covarianceOf(const TObject* h) {
  auto c=tag(h,"physics_covariance");
  return c.is_null()?unknown("upstream covariance availability was not recorded"):c;
}
inline Trace derivedTrace(const std::string& operation,const std::vector<const TObject*>& inputs,const Json& parameters=Json::object()) {
  std::vector<Trace> traces;
  for(const auto* object:inputs)traces.push_back(object?projectionOf(object,Json::array()):Trace::unavailable("missing operand"));
  return Trace::combine(operation,traces,parameters);
}
inline Json mapValues(const TH1* h,bool errorsComputed) {
  Json values=histogram(h);
  values["errors"]=errorsComputed?"native_variances":"not_computed";
  if(!errorsComputed)for(auto& cell:values["cells"])cell["variance"]=nullptr;
  return values;
}
inline Trace selectedTrace(const TH1* h,double low,double high) {
  const Json window=Json::array({{{"axis",0},{"low",low},{"high",high}}});
  return projectionOf(h,window).then("select_complete_bins",{{"low",low},{"high",high}});
}

/* Bind the numerical object before style-only marker offsets or line sampling
   can become a scientific value. Snapshot checking catches a stale binding
   after an unrecorded SetPoint/SetBinContent mutation. */
inline std::string seriesRole(const TObject* object,const std::string& fallback) {
  const auto named=tag(object,"physics_series_role");
  if(named.is_string())return named.get<std::string>();
  const auto binding=tag(object,"physics_series");
  return binding.is_object()?binding.at("role").get<std::string>():fallback;
}
inline void nameSeries(TObject* object,const std::string& role) {
  require(!role.empty(),"series slot needs a name");setTag(object,"physics_series_role",role);
}

template<class T> T* bind(T* object,const std::string& role,const Trace& trace,const Json& value,
                       const Json& covarianceRecord) {
  require(object && !role.empty(),"series role required");
  attachTrace(object,trace);setTag(object,"physics_covariance",covarianceRecord);
  setTag(object,"physics_series",{{"role",role},{"trace",trace.json()},{"value",portable(value)},
    {"covariance",covarianceRecord},{"drawing_snapshot",nativeValues(object)}});
  return object;
}
template<class T> T* bindNative(T* object,const std::string& role,const Trace& trace,
                             const Json& cov=unknown("producer did not declare covariance availability")) {
  return photonjet_science::bind(object,role,trace,nativeValues(object),cov);
}
template<class T> T* bindMap(T* object,const std::string& role,const Trace& trace,bool errorsComputed) {
  attachTrace(object,trace);
  return photonjet_science::bind(object,role,trace,mapValues(object,errorsComputed),
      covariance("not_computed",{},0,"native map calculation does not propagate cross-cell covariance"));
}
inline void bindFit(TF1* fit,const std::string& role,const Trace& input,const TFitResultPtr& result,const Json& method) {
  Json cov=unknown("fit did not return a covariance result");
  if(result.Get()) {
    const int n=result->NPar();std::vector<double> values;
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)values.push_back(result->CovMatrix(i,j));
    cov=covariance("measured",values,n,"ROOT fit covariance; status "+std::to_string(int(result)));
  }
  auto value=nativeValues(fit);value["fit_status"]=int(result);
  value["covariance_status"]=result.Get()?result->CovMatrixStatus():-1;
  if(!result.Get() || int(result)!=0)value["state"]="missing";
  photonjet_science::bind(fit,role,input.then("ROOT_fit",method),value,cov);
}
inline void decoration(TObject* object,const std::string& reason) {
  require(!reason.empty(),"decoration exclusion needs a reason");
  setTag(object,"physics_decoration",{{"reason",reason}});
}
inline void copyBinding(TObject* out,const TObject* input) {
  const auto b=tag(input,"physics_series");
  if(!b.is_null()) { auto copy=b;copy["drawing_snapshot"]=nativeValues(out);setTag(out,"physics_series",copy); }
}
inline Json scalar(double value,const std::string& unit,const std::string& method,
                   const Json& uncertainty=absent("no uncertainty calculated for this annotation")) {
  require(!unit.empty() && !method.empty(),"scalar needs units and method");
  return {{"type","scalar"},{"state",std::isfinite(value)?"measured":"missing"},
    {"value",std::isfinite(value)?Json(real(value)):Json()},{"unit",unit},{"method",method},{"uncertainty",uncertainty}};
}
inline std::string bindingKey(const Json& binding);
inline void expect(TVirtualPad& pad,const Json& binding);
inline void annotation(TVirtualPad& pad,const std::string& role,const Json& value,const Trace& trace) {
  require(!role.empty(),"annotation role required");
  const Json binding={{"role",role},{"trace",trace.json()},{"value",portable(value)},
    {"covariance",absent("scalar uncertainty is recorded in the annotation")}};
  expect(pad,binding);
  pad.GetListOfPrimitives()->Add(new TNamed("physics_annotation",binding.dump().c_str()));
}
inline void attachAnnotation(TObject* object,const std::string& role,const Json& value,const Trace& trace) {
  require(!role.empty(),"annotation role required");
  Json rows=tag(object,"physics_annotations");if(!rows.is_array())rows=Json::array();
  rows.push_back({{"role",role},{"trace",trace.json()},{"value",portable(value)},
    {"covariance",absent("scalar uncertainty is recorded in the annotation")}});
  setTag(object,"physics_annotations",rows);
}
inline std::string bindingKey(const Json& binding) {
  return Json({{"role",binding.at("role")},{"channels",Trace::fromJson(binding.at("trace")).identity()}}).dump();
}
inline void expect(TVirtualPad& pad,const Json& binding) {
  pad.GetListOfPrimitives()->Add(new TNamed("physics_expected",bindingKey(binding).c_str()));
}
/* Keep the output's declared membership independently of its primitives. A
   subsequent accidental removal/clear can no longer produce a complete but
   silently shortened receipt. Frames/legend-only glyphs are excluded only by
   their explicit decoration tag, never by an object name or pixel property. */
template<class T> void draw(T* object,const char* option="") {
  require(object,"null object passed to native draw");
  auto* pad=gPad;
  std::vector<TObject*> declarations;
  if(pad)for(auto* item:*pad->GetListOfPrimitives())
    if(std::string(item->GetName())=="physics_annotation" || std::string(item->GetName())=="physics_expected")declarations.push_back(item);
  if(pad)for(auto* item:declarations)pad->GetListOfPrimitives()->Remove(item);
  object->Draw(option);
  if(pad)for(auto* item:declarations)pad->GetListOfPrimitives()->Add(item);
  if(pad && (dynamic_cast<TH1*>(object)||dynamic_cast<TGraph*>(object)||dynamic_cast<TF1*>(object)) &&
      tag(object,"physics_decoration").is_null()) {
    const auto binding=tag(object,"physics_series");
    if(!binding.is_null())expect(*pad,binding);
    const auto annotations=tag(object,"physics_annotations");
    if(annotations.is_array())for(const auto& value:annotations)expect(*pad,value);
  }
}
/* Some native renderers paint a histogram as individual TBoxes. Retain its
   scientific binding on the pad without drawing an additional histogram. */
inline void observe(TVirtualPad& pad,const TObject* object) {
  const auto binding=tag(object,"physics_series");
  if(binding.is_null()) {
    annotation(pad,"unbound_map",unknown("manual painter received an unbound scientific object"),traceOf(object));return;
  }
  require(binding.at("drawing_snapshot")==nativeValues(object),"manual painter received stale numerical binding");
  expect(pad,binding);
  pad.GetListOfPrimitives()->Add(new TNamed("physics_annotation",binding.dump().c_str()));
}

/* ROOT DrawFrame clears the pad. Typed annotations are declarations of the
   output's measurements and must survive that layout operation. Preserve only
   these explicitly typed records, never arbitrary old drawn primitives. */
inline TH1F* drawFrame(TVirtualPad& pad,double xlow,double ylow,double xhigh,double yhigh) {
  std::vector<TObject*> declarations;
  for(TObject* object:*pad.GetListOfPrimitives())
    if(std::string(object->GetName())=="physics_annotation" || std::string(object->GetName())=="physics_expected")declarations.push_back(object);
  for(auto* object:declarations)pad.GetListOfPrimitives()->Remove(object);
  auto* frame=pad.DrawFrame(xlow,ylow,xhigh,yhigh);
  decoration(frame,"ROOT axis frame created without data");
  for(auto* object:declarations)pad.GetListOfPrimitives()->Add(object);
  return frame;
}

/* Fresh collector per output. The archive/provenance lives outside the
   PlotPhysics identity. A family can only be admitted when both the native
   fixture manifest and this output's complete measurement roster pass. */
struct Collector {
  Json bindings=Json::object(), missing=Json::array(), origins=Json::object();
  std::set<const TObject*> seen;
  std::set<std::string> expected;
  void add(const Json& b) {
    const Trace t=Trace::fromJson(b.at("trace"));
    const std::string key=bindingKey(b);  // no values, windows or definition content
    Json scientific=b;scientific.erase("drawing_snapshot");scientific["trace"]=t.json();scientific["trace"].erase("origins");
    for(const auto& [id,origin]:t.origins.items()) { require(!origins.contains(id)||origins.at(id)==origin,"conflicting canonical origins");origins[id]=origin; }
    if(bindings.contains(key) && bindings.at(key)!=scientific)missing.push_back("conflicting values for stable series role: "+b.at("role").get<std::string>());
    else bindings[key]=scientific;
  }
  void visit(TPad& pad) {
    for(TObject* object:*pad.GetListOfPrimitives()) {
      if(!seen.insert(object).second)continue;
      if(auto* child=dynamic_cast<TPad*>(object)){visit(*child);continue;}
      if(std::string(object->GetName())=="physics_expected") { expected.insert(object->GetTitle());continue; }
      if(std::string(object->GetName())=="physics_annotation") { add(Json::parse(object->GetTitle()));continue; }
      if(!dynamic_cast<TH1*>(object) && !dynamic_cast<TGraph*>(object) && !dynamic_cast<TF1*>(object))continue;
      if(!tag(object,"physics_decoration").is_null())continue;
      const auto annotations=tag(object,"physics_annotations");
      if(annotations.is_array())for(const auto& value:annotations)add(value);
      const auto b=tag(object,"physics_series");
      if(b.is_null()) { missing.push_back(std::string("unbound scientific primitive: ")+object->ClassName());continue; }
      if(b.at("drawing_snapshot")!=nativeValues(object))missing.push_back("scientific object changed after its binding");
      add(b);
    }
  }
  Json finish(TPad& pad,const std::string& output,const std::string& family,bool annotationRosterComplete) {
    visit(pad);
    for(const auto& key:expected)if(!bindings.contains(key))missing.push_back("declared scientific series missing: "+key);
    photonjet_render::PlotPhysics result(family);
    std::map<std::string,bool> complete;
    for(const auto& c:categories())complete[c]=!bindings.empty() && missing.empty();
    bool quantities=!bindings.empty() && missing.empty(),covarianceComplete=quantities;
    for(const auto& [key,b]:bindings.items()) {
      const Trace t=Trace::fromJson(b.at("trace"));
      result.quantity(key,b.at("value"));result.definition(key,t.definitions());
      result.dependency(key,t.scientificJson().at("roots"));
      result.definition(key+"/covariance",b.at("covariance"));
      for(const auto& c:categories())complete[c]=complete[c] && t.complete(c);
      quantities=quantities && measurementsComplete(b.at("value"));
      const std::string state=b.at("covariance").value("state","");
      covarianceComplete=covarianceComplete && (state=="measured" || state=="not_applicable" || state=="not_computed");
    }
    if(bindings.empty())result.definition("missing",{{"reason","no bound scientific values"}});
    for(const auto& c:categories())result.coverage(c,complete[c]?"complete":"unknown","typed per-series producer trace and operations");
    result.coverage("quantities",quantities?"complete":"unknown","exact values and explicit measurement availability");
    result.coverage("errors",quantities?"complete":"unknown","recorded native variances or asymmetric error arrays");
    result.coverage("covariance",covarianceComplete?"complete":"unknown","actual covariance or producer-declared calculation availability");
    result.coverage("scalar_annotations",annotationRosterComplete && missing.empty()?"complete":"unknown","family annotation roster and explicit typed scalars");
    bool dependencies=!bindings.empty() && missing.empty();for(const auto& [c,value]:complete)dependencies=dependencies&&value;
    result.coverage("dependencies",dependencies?"complete":"unknown","only selected typed per-output dependencies");
    Json failures=Json::array();
    for(const auto& [key,b]:bindings.items()) {
      const auto t=Trace::fromJson(b.at("trace"));Json unknownCategories=Json::array();
      for(const auto& category:categories())if(!t.complete(category))unknownCategories.push_back(category);
      if(!unknownCategories.empty()||!measurementsComplete(b.at("value")))
        failures.push_back({{"role",b.at("role")},{"unknown_categories",unknownCategories},{"missing",t.missing},{"measurements_complete",measurementsComplete(b.at("value"))}});
    }
    Json receipt=result.finish(output);receipt["diagnostics"]={{"unbound",missing},{"bound_series",bindings.size()},{"incomplete",failures}};return receipt;
  }
};
}  // namespace photonjet_science
#endif
