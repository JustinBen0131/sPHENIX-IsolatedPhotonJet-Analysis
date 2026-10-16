#pragma once
/* Source spans are provenance, never scientific identity. The native settings
   reader supplies each evaluated value; no Python parser guesses C++ types. */
#include "RenderMetadata.h"
#include "../../TreeToHists/src/ScientificTrace.h"
#include <yaml.h>
#include <filesystem>
#include <set>
#include <cstdlib>
#include <fstream>
#include <regex>
#include <sstream>

namespace fa {
struct ConfigurationScientificWitness {
  using Json=photonjet_science::Json;
  Json sites=Json::object(),exports=Json::object();
  std::string source,relative,sourceSha;
  std::set<std::string> ambiguous;

  static std::string relativeSource(const std::filesystem::path &path) {
    /* Keep the actual input location. A run configuration outside the admitted
       source tree cannot impersonate canonical config/analysis.yaml. */
    auto absolute=std::filesystem::absolute(path).lexically_normal();
    if (const char *root=std::getenv("PHOTONJET_SCIENTIFIC_SOURCE_ROOT")) {
      auto rel=absolute.lexically_relative(std::filesystem::absolute(root).lexically_normal());
      if(!rel.empty() && *rel.begin()!="..")return rel.generic_string();
    }
    for(auto p=absolute.parent_path();!p.empty() && p!=p.parent_path();p=p.parent_path())
      if(p.filename()=="FinalAnalysis")return absolute.lexically_relative(p.parent_path()).generic_string();
    return absolute.generic_string();
  }
  void span(const std::string &key,size_t first,size_t last) {
    if(first>=last || last>source.size())return;
    if(sites.contains(key) || ambiguous.count(key)) {sites.erase(key);ambiguous.insert(key);return;}
    size_t radius=48;
    std::string before,after;
    for(;;) {
      const auto begin=first>radius?first-radius:0,finish=std::min(source.size(),last+radius);
      before=source.substr(begin,first-begin);after=source.substr(last,finish-last);
      const auto exact=source.substr(begin,finish-begin);
      if(source.find(exact)==source.rfind(exact))break;
      radius*=2;
    }
    if(before.empty() || after.empty())return;
    sites[key]={{"kind","source_text"},{"path",relative},{"source_sha256",sourceSha},{"encoding","native_value"},
      {"canonical_key",key},{"before_text",source.substr(first,last-first)},{"selector",{{"before_context",before},{"after_context",after}}}};
  }
  void yamlNode(yaml_document_t &document,yaml_node_t *node,const std::string &prefix) {
    if(!node)return;
    if(node->type==YAML_MAPPING_NODE) {
      for(auto pair=node->data.mapping.pairs.start;pair<node->data.mapping.pairs.top;++pair) {
        auto *key=yaml_document_get_node(&document,pair->key),*value=yaml_document_get_node(&document,pair->value);
        if(!key || key->type!=YAML_SCALAR_NODE)continue;
        const std::string text(reinterpret_cast<const char*>(key->data.scalar.value),key->data.scalar.length);
        if(text=="<<")continue; // merged/aliased definitions require explicit admission, never a guessed edit site
        const std::string name=prefix.empty()?text:prefix+"."+text;
        if(value && value->start_mark.index<key->end_mark.index)continue;
        yamlNode(document,value,name);
      }
    } else if(!prefix.empty()) span(prefix,node->start_mark.index,node->end_mark.index);
  }
  void load(const std::string &path) {
    std::ifstream file(path,std::ios::binary);
    source.assign(std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>());
    relative=relativeSource(path);sourceSha=photonjet_render::sha256(source);
    if(source.empty())return;
    const auto extension=std::filesystem::path(path).extension();
    if(extension==".yaml" || extension==".yml") {
      yaml_parser_t parser{};yaml_document_t document{};
      if(!yaml_parser_initialize(&parser))return;
      yaml_parser_set_input_string(&parser,reinterpret_cast<const unsigned char*>(source.data()),source.size());
      if(yaml_parser_load(&parser,&document)) {yamlNode(document,yaml_document_get_root_node(&document),"");yaml_document_delete(&document);}
      yaml_parser_delete(&parser);
    } else {
      const std::regex line(R"((^|\n)([A-Za-z][A-Za-z0-9_.]*)[ \t]*:[ \t]*([^#\r\n]+))");
      for(auto i=std::sregex_iterator(source.begin(),source.end(),line);i!=std::sregex_iterator();++i) {
        const auto &match=*i;const std::string value=match[3];
        const size_t length=value.find_last_not_of(" \t")+1;
        span(match[2],match.position(3),match.position(3)+length);
      }
    }
  }
  Json evaluated(const std::string &key,const Json &value) {
    if(!sites.contains(key))return Json();
    Json site=sites.at(key);
    const std::string id="FinalAnalysis.setting/"+key+"/"+sourceSha;
    site["export_id"]=id;
    exports[id]={{"path",relative},{"source_sha256",sourceSha},{"selector",site.at("selector")},
      {"canonical_key",key},{"before_text",site.at("before_text")},{"value",photonjet_science::portable(value)}};
    return site;
  }
};
}  // namespace fa
