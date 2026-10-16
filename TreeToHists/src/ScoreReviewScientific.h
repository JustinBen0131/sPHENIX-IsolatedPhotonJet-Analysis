/* Scientific metadata for the existing ScoreReview fill/collection. */
#pragma once
#include "QAScientificContext.h"
#include "../../FinalAnalysis/src/RootScientificTrace.h"
namespace photonjet_science {
inline Json scoreReviewDefinition(const Json& context,int centralitySlot,const std::string& population) {
  Json terms=Json::array(),declarations=Json::object();
  for(const auto& category:categories())declarations[category]=known("native score-review term roster");
  const auto add=[&](const std::string& category,const std::string& key,const Json& value){terms.push_back(term(category,key,value));};
  for(const auto* key:{"collision","event_record_selection","ownership","candidate_geometry","preselection","truth_classification"})
    add("selection",key,child(context,"/"+std::string(key)));
  add("selection","review_population",{{"et",{15.,35.}},{"eta_abs_lt",.7},{"vertex_abs_lt_cm",10.},
    {"centrality",centralitySlot==0?Json({0.,20.}):centralitySlot==1?Json({20.,50.}):Json({50.,80.})},
    {"finite_reconstructed_SUB1_R04_witness",true},{"isolation_cut",false},
    {"signal_winner","matched source-role-1 truth; smallest metric, largest ET, lowest index"},
    {"background","all eligible inclusive candidates; truth_background uses unambiguous nonprompt dominance"}});
  add("selection","isolation_witness",child(context,"/isolation_observable"));
  add("models","npb",child(context,"/ncb_model"));add("models","score",child(context,"/score_model"));
  add("weights","event",child(context,"/base_event_weight"));
  add("normalization","sample_factor",{{"applied",false}});
  add("projection","score_accumulation",{{"bins",50},{"range",{0.,1.}},
      {"exact_one","nextafter(1,0)"},{"content","sum(w_event)"},{"variance","sum(w_event*w_event)"}});
  return definition("score_review",{{"centrality_slot",centralitySlot},{"population",population},
      {"source_channel",context.value("source_channel",Json::object())}},terms,declarations);
}
inline Trace scoreReviewNormalization(const Trace& input,const std::string& sample,double crossSection,long long denominator) {
  Json declarations=Json::object();for(const auto& category:categories())declarations[category]=known("native sample normalization");
  const auto normalization=Trace::source(definition("score_review_sample_normalization",{{"sample",sample}},
    Json::array({term("normalization","factor",{{"cross_section_pb",crossSection},{"denominator",denominator},
      {"factor",crossSection/denominator},{"applied_once",true}})}),declarations));
  return Trace::combine("apply_sample_normalization",{input,normalization});
}
}
