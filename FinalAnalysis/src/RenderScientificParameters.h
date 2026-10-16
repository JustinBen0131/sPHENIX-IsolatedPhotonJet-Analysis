#pragma once
#include "RootScientificTrace.h"
namespace photonjet_science {
/* Exact settings read by the producer, compared to the value actually used by
   the renderer. Missing or inconsistent witnesses cannot certify a drawing. */
inline Trace renderSetting(TDirectory& file,const std::string& key,const Json& effective) {
  auto* object=dynamic_cast<TObjString*>(file.Get("scientific_parameter_bindings"));
  if(!object)object=dynamic_cast<TObjString*>(file.Get("meta/scientific_parameter_bindings")); // early capture compatibility
  if(!object)return Trace::unavailable("native setting witness absent: "+key);
  const auto records=Json::parse(object->GetString().Data());
  const auto global=records.value("global",Json::object());
  if(!global.contains(key))return Trace::unavailable("native setting key absent: "+key);
  const auto& binding=global.at(key);
  if(portable(binding.at("value"))!=portable(effective))return Trace::unavailable("native setting disagrees with effective render value: "+key);
  return Trace().withParameter(binding.at("reference"),binding.at("category"),effective,binding.value("origin",Json()));
}
inline void renderSetting(TVirtualPad& pad,TDirectory& file,const std::string& key,const Json& effective) {
  annotation(pad,"renderer_setting/"+key,{{"type","effective_render_setting"},{"value",portable(effective)}},renderSetting(file,key,effective));
}
inline Trace savedScalarTrace(TDirectory& file,const std::string& path) {
  auto* object=dynamic_cast<TObjString*>(file.Get(path.c_str()));
  if(!object)return Trace::unavailable("native scalar witness absent: "+path);
  const auto record=Json::parse(object->GetString().Data());
  return record.contains("trace")?Trace::fromJson(record.at("trace")):Trace::unavailable("native scalar trace absent: "+path);
}
} // namespace photonjet_science
