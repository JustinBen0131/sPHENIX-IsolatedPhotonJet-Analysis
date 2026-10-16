/* Public source-bound cut plans and grouped reduction. Included by TreeToHists.cc
   so all normal plan, inventory, source seal and merge functions remain shared.
   No scheduler, calibration derivation or alternate physics reducer lives here. */
namespace treetohists {
namespace selection_fanout {
using photonjet::selection_views::require;
// A grid's ordinary reference is historical nominal; its cut-view products
// must instead carry the explicitly selected corrected nominal response.
std::string viewRecipe(const Json& parent){
 if(!parent.contains("response_component_grid"))return parent.at("histograms").at("response_variation");
 for(const auto& sample:parent.at("samples"))require(sample.at("family")!="data","grid selection plans require a SIM-only parent");
 for(const auto& variant:parent.at("response_component_grid").at("variants"))if(variant.at("name")=="smnom"){
  const std::string text=variant.at("variation");const auto parsed=photonjet::response::parseVariation(text.c_str());
  require(photonjet::response::responseOnly(parsed) && parsed.photonResolution=="ppg18_nominal" &&
   parsed.jetResolution=="sphenix_pp2026_nominal" && parsed.seed==480 && !parsed.centralityOff &&
   parsed.jetScale==1. && parsed.emScale==1.,"smnom grid component differs from the selected nominal recipe");
  return text;
 }
 throw std::runtime_error("selection grid lacks its smnom reference component");
}
Json sealed(Json value,const char* key="plan_sha256"){value.erase(key);value[key]=digest(value);return value;}
Json checkedPlan(const fs::path& work){auto p=readJson((work/"plan.json").string());auto b=p;b.erase("plan_sha256");require(digest(b)==p.at("plan_sha256").get<std::string>(),"plan checksum differs");checkReductionCode(p);return p;}
Json groupPlan(const fs::path& work,const std::string& group){
 const auto shard=std::to_string(std::stoul(group.substr(1))/1000);
 auto d=readJson((work/"fill_groups"/shard/(group+".json")).string());auto body=d;body.erase("descriptor_sha256");
 require(digest(body)==d.at("descriptor_sha256").get<std::string>() && d.at("plan").at("plan_sha256")==readJson((work/"binding.json").string()).at("plan_sha256"),"group binding differs");return d;
}
void localGroups(const fs::path& work){
 auto p=checkedPlan(work);require(!fs::exists(work/"fill_groups")&&!fs::exists(work/"merge_rows"),"group graph already exists");
 size_t groupId=0,nodeId=0;Json graph=Json::array();
 for(const auto& sample:p.at("samples")){
  Json jobs=Json::array();for(const auto& j:p.at("jobs"))if(j.at("sample")==sample.at("key"))jobs.push_back(j);
  std::vector<std::string> inputs;Json inputSamples=Json::array();
  if(!jobs.empty()){
   const std::string id="G"+std::to_string(groupId++),shard=std::to_string((groupId-1)/1000);auto sub=sample;sub["intermediate"]=true;
   Json small={{"schema",p.at("schema")},{"producer_code",p.at("producer_code")},{"selection",p.at("selection")},{"histograms",p.at("histograms")},{"root_version",p.at("root_version")},{"plan_sha256",p.at("plan_sha256")},{"samples",Json::array({sub})},{"jobs",jobs}};
   if(p.contains("response_component_grid"))small["response_component_grid"]=p.at("response_component_grid");
   Json d=sealed(Json{{"plan",small}},"descriptor_sha256");fs::create_directories(work/"fill_groups"/shard);writeText(work/"fill_groups"/shard/(id+".json"),d.dump()+"\n");
   inputs.push_back((work/"chunks"/shard/id/(sample.at("key").get<std::string>()+".root")).string());inputSamples.push_back(sub);graph.push_back({{"id",id},{"operation","fill-selection-group"}});
  }else if(sample.contains("run_packages"))for(const auto& key:sample.at("run_packages")){
   inputs.push_back((work/"final"/(key.get<std::string>()+".root")).string());inputSamples.push_back(sampleOf(p,key));
  }else throw std::runtime_error("selection groups: source jobs absent");
  const std::string id="M"+std::to_string(nodeId++);Json small=p;small["jobs"]=Json::array();
  Json node={{"sample",sample},{"plan",small},{"inputs",inputs},{"input_samples",inputSamples},{"output_directory",(work/"final").string()},{"plan_sha256",p.at("plan_sha256")}};node=sealed(node,"descriptor_sha256");fs::create_directories(work/"merge_rows");writeText(work/"merge_rows"/(id+".json"),node.dump()+"\n");graph.push_back({{"id",id},{"operation","merge-group"}});
 }
 writeText(work/"selection_graph.json",Json{{"schema","PhotonJetNativeSelectionGraphV2"},{"plan_sha256",p.at("plan_sha256")},{"nodes",graph},{"submission_authorized",false}}.dump(2)+"\n");
}
void planView(const fs::path& parentWork,const fs::path& output,const std::string& prescription,const std::string& source){
 require(!fs::exists(output),"view work already exists; reuse or reconcile it");auto p=checkedPlan(parentWork);
 require(!p.contains("selection_fanout"),"selection parent must be an ordinary complete reference plan");
 require(!p.at("histograms").value("candidate_records",false),"native cut fanout currently requires the explicitly selected no-candidate-archive setting; event bootstrap records are retained");
 const auto selectedRecipe=viewRecipe(p);
 const auto recipe=photonjet::response::parseVariation(selectedRecipe.c_str());
 require(photonjet::selection_views::cutComposition(recipe),"only nominal or the exact corrected nominal resolution recipe can accompany cuts");
 const bool jes=prescription=="data_jes_up" || prescription=="data_jes_down";
 Json contract;
 if(jes){
  require(p.at("histograms").at("data_jet_scale").at("pp")==1./.916 && p.at("histograms").at("data_jet_scale").at("auau")==1.,"JES parent must use exact nominal scales");
  contract={{"schema","PhotonJetNativeSelectionPlanV2"},{"kind","DATA_JES"},{"prescription",prescription},{"parent_selection_sha256",p.at("selection").at("selection_sha256")},{"effective_selection_sha256",p.at("selection").at("selection_sha256")},{"selection",p.at("selection")},{"definition_sources",Json::array()},{"admission","EXISTING_NATIVE_DATA_JES_FUSION"}};
  p["histograms"]["data_jet_scale"]["pp"]=prescription=="data_jes_up"?1./.907:1./.923;
 }else contract=photonjet::selection_views::adopted(p.at("selection"),prescription,source,digest,sha256File);
 p["histograms"]["response_variation"]=selectedRecipe;
 contract["parent_plan_sha256"]=p.at("plan_sha256");contract["reference_response_variation"]=selectedRecipe;
 p.erase("response_component_grid");p["selection"]=contract.at("selection");contract.erase("selection");p["selection_fanout"]=contract;p=sealed(p);savePlan(p,output.string());
 if(fs::exists(parentWork/"fill_groups"))for(const auto& entry:fs::recursive_directory_iterator(parentWork/"fill_groups"))if(entry.is_regular_file()&&entry.path().extension()==".json"){
  auto d=readJson(entry.path().string());auto& group=d["plan"];group.erase("response_component_grid");group["selection"]=p.at("selection");group["histograms"]=p.at("histograms");group["selection_fanout"]=contract;group["plan_sha256"]=p.at("plan_sha256");d=sealed(d,"descriptor_sha256");
  const auto target=output/"fill_groups"/fs::relative(entry.path(),parentWork/"fill_groups");fs::create_directories(target.parent_path());writeText(target,d.dump()+"\n");
 }
 if(fs::exists(parentWork/"merge_rows"))for(const auto& entry:fs::directory_iterator(parentWork/"merge_rows"))if(entry.is_regular_file()&&entry.path().extension()==".json"){
  auto node=readJson(entry.path().string());node["plan_sha256"]=p.at("plan_sha256");node["plan"]["plan_sha256"]=p.at("plan_sha256");node["plan"].erase("response_component_grid");node["plan"]["selection"]=p.at("selection");node["plan"]["histograms"]=p.at("histograms");node["plan"]["selection_fanout"]=contract;
  const auto remap=[&](const std::string& path){const auto from=parentWork.string();return path.rfind(from+"/",0)==0?output.string()+path.substr(from.size()):path;};
  node["output_directory"]=remap(node.at("output_directory"));for(auto& input:node["inputs"])input=remap(input.get<std::string>());
  node=sealed(node,"descriptor_sha256");fs::create_directories(output/"merge_rows");writeText(output/"merge_rows"/entry.path().filename(),node.dump()+"\n");
 }
 writeText(output/"selection_fanout.json",contract.dump(2)+"\n");
}
struct Products {
 photonjet::Profile profile;photonjet::Histograms measurement,timing;photonjet::EventTable records,timed;
 photonjet::qa::Book qa;photonjet::response::ResponseAccumulators response;
 Products(photonjet::Profile p,const std::vector<photonjet::CentralityClass>& classes,const photonjet::response::ResponseGrids& grids,bool sim,photonjet::qa::Family family)
  :profile(std::move(p)),measurement(classes),timing(classes),qa(true){
   photonjet::qa::prebook(qa,family,profile.system,classes);
   if(sim){response.resize(classes.size());for(size_t c=0;c<classes.size();++c)for(double t:photonjet::kJetThresholds)response[c].emplace_back("_"+classes[c].name+"_"+photonjet::thresholdName(t),grids,profile.system==photonjet::System::pp);}
 }
};
void configure(const Json& p,const Json& sample,const Json& job){
 checkReductionCode(p);checkVertexFit(sample);const auto periods=ppPeriodsOf(p,sample);checkPpPeriodFits(periods);
 require(p.at("root_version")==gROOT->GetVersion(),"ROOT runtime differs");const auto& h=p.at("histograms");
 photonjet::configuredXjEdges=h.at("native_xj_edges").get<std::vector<double>>();photonjet::configuredMeasurementEtEdges=h.at("native_measurement_et_edges_gev").get<std::vector<double>>();
 photonjet::configureJetRadius(h.at("jet_radius"));photonjet::configuredDataJetScale=h.at("data_jet_scale").at(sample.at("system")=="PP"?"pp":"auau");
 photonjet::candidate_cache::enabled=false;photonjet::strictReduction=true;photonjet::configurePpPeriods(periods.is_null()?Json::object():periods);
 photonjet::configureVertexFit(sample.at("vertex_reweighting").is_null()?Json::object():sample.at("vertex_reweighting"));
 photonjet::auauReceiptSha=sample.at("vertex_reweighting").is_null()?"":sample.at("vertex_reweighting").at("evidence_sha256").get<std::string>();photonjet::auauStitchingSha=sample.at("normalization").value("stitching_table_sha256","");
 photonjet::configuredOwnership=ownershipOf(sample.at("normalization"));
 const auto& source=sample.at("sources").at(job.at("source").get<size_t>());photonjet::entryBegin=job.value("range_begin",0LL);photonjet::entryEnd=job.value("range_end",source.at("metadata").at("entries").get<long long>());
 photonjet::progressRowTotal=photonjet::entryEnd-photonjet::entryBegin;photonjet::progressRowProcessed=0;photonjet::progressOffset=0;
 photonjet::progressTotal=photonjet::progressRowTotal*(sample.at("normalization").contains("ownership")?4:3);photonjet::progressStart=std::chrono::steady_clock::now();
}
void row(const std::vector<Json>& plans,const std::vector<fs::path>& works,size_t index,double mbdReference){
 const auto& parent=plans[0];const auto& job=parent.at("jobs").at(index);const auto& sample=parent.at("samples").at(0);const auto& source=sample.at("sources").at(job.at("source").get<size_t>());
 configure(parent,sample,job);const std::string input=job.at("input");const auto stamp=std::make_pair(fs::file_size(input),fs::last_write_time(input));
 std::string familyName,systemName;const auto observed=sourceInfo(input,parent.at("selection"),familyName,systemName);require(observed.at("identity")==source.at("identity")&&observed.at("metadata")==source.at("metadata"),"source changed since planning");
 const auto system=sample.at("system")=="PP"?photonjet::System::pp:photonjet::System::auau;const bool sim=sample.at("family")!="data";
 const auto family=!sim?photonjet::qa::Family::data:sample.at("family")=="photonjet"?photonjet::qa::Family::phosim:photonjet::qa::Family::jetsim;
 const auto classes=photonjet::centralityClasses(system);const auto& h=parent.at("histograms");photonjet::response::ResponseGrids grids;
 grids.reco=h.at("native_reco_response_edges_gev").get<std::vector<double>>();grids.truth=h.at("native_truth_response_edges_gev").get<std::vector<double>>();grids.xj=h.at("native_xj_edges").get<std::vector<double>>();
 const auto weights=sim&&system==photonjet::System::auau?photonjet::centralityWeights(listed(sample.at("weights")),system):photonjet::CentralityWeights{};
 const auto recipe=photonjet::response::parseVariation(h.at("response_variation").get<std::string>().c_str());require(photonjet::selection_views::cutComposition(recipe),"unsupported cut composition");
 std::vector<std::unique_ptr<Products>> products;std::vector<fs::path> dirs;std::vector<photonjet::response::Variation> recipes;
 for(size_t v=0;v<plans.size();++v){
  recipes.push_back(photonjet::response::parseVariation(plans[v].at("histograms").at("response_variation").get<std::string>().c_str()));
  dirs.emplace_back(jobDirectory(works[v].string(),plans[v].at("jobs").at(index)));fs::create_directories(dirs.back());
  auto profile=photonjet::configuredProfile("nominal",system,(works[v]/"selection.json").c_str());products.push_back(std::make_unique<Products>(profile,classes,grids,sim,family));
 }
 const auto isJes=[&](size_t v){return v && plans[v].at("selection_fanout").value("kind","")=="DATA_JES";};
 std::map<size_t,std::unique_ptr<photonjet::DataJesView>> jesProducts;
 for(size_t v=1;v<plans.size();++v)if(isJes(v)){require(!sim && system==photonjet::System::pp,"JES views require pp DATA");jesProducts.emplace(v,std::make_unique<photonjet::DataJesView>(plans[v].at("histograms").at("data_jet_scale").at("pp"),classes));}
 const auto expectedGap=[&](size_t v){return v&&!isJes(v)&&system==photonjet::System::pp?plans[v].at("selection_fanout").value("expected_parent_stored_gap_gev",std::numeric_limits<double>::quiet_NaN()):std::numeric_limits<double>::quiet_NaN();};
 const std::string sourceText="["+photonjet::sourceRecord(input)+"]";
 Json gridReceipt;
 photonjet::progressPhase="selection_measurement_response";
 if(!sim){
  std::vector<photonjet::SelectionMeasurementView> values;std::vector<photonjet::SelectionMeasurementView*> views;values.reserve(products.size()-1);
  for(size_t v=1;v<products.size();++v)if(!isJes(v)){photonjet::SelectionMeasurementView view;view.profile=products[v]->profile;view.histograms=&products[v]->measurement;view.records=&products[v]->records;view.timing=&products[v]->timing;view.timingRecords=&products[v]->timed;view.expectedStoredGap=expectedGap(v);values.push_back(view);}
  for(auto& value:values)views.push_back(&value);
  std::vector<photonjet::DataJesView*> sharedJes;for(auto& item:jesProducts)sharedJes.push_back(item.second.get());
  photonjet::fillFile(input,0,products[0]->profile,classes,products[0]->measurement,products[0]->records,{},"",0,nullptr,&products[0]->timing,&products[0]->timed,nullptr,"",{}, sharedJes,views);
 }else{
  std::vector<photonjet::response::ResponseView> views;
  for(size_t v=0;v<products.size();++v){photonjet::response::ResponseView view;view.variation=recipes[v];view.acc=&products[v]->response;view.selection=v?&products[v]->profile:nullptr;view.measurement=&products[v]->measurement;view.measurementRecords=&products[v]->records;view.expectedStoredGap=expectedGap(v);views.push_back(view);}
  if(parent.contains("response_component_grid")){
   namespace ng=photonjet::response::native_grid;const auto& config=parent.at("response_component_grid");ng::validateSettings(config);
   Json spec={{"schema","PhotonJetNativeResponseGridV1"},{"identity",{{"task","THE-510"},{"workstream","selection-grid"},{"campaign",parent.at("plan_sha256")},{"row",std::to_string(index)}}},
    {"sources",Json::array({{{"path",input},{"sha256",source.at("expected_sha256")},{"input_file_sha256",source.at("identity")},{"entries",source.at("metadata").at("entries")},{"entry_begin",photonjet::entryBegin},{"entry_end",photonjet::entryEnd}}})},
    {"variants",config.at("variants")},{"producer_code",parent.at("producer_code")},{"runtime_context",ng::context(products[0]->profile,grids,weights)},
    {"limits",config.at("limits")},{"output_directory",fs::absolute(dirs[0]/"response_grid").lexically_normal().string()}};
   ng::copyEmsContract(config,spec);
   for(auto& v:spec["variants"]){v["plan"]=ng::resolvedPlan(spec,v);v["plan_sha256"]=digest(v["plan"]);}
   gridReceipt=ng::run(spec,products[0]->profile,grids,weights,views);
   require(gridReceipt.at("original_entries")==photonjet::progressRowTotal,"grid ended before full row");
   for(const auto& product:gridReceipt.at("products")){
    const std::string name=product.at("name");const auto target=dirs[0]/(name=="nominal"?"response.root":"response__"+name+".root");
    fs::create_hard_link(dirs[0]/"response_grid"/product.at("path").get<std::string>(),target);
   }
  }else photonjet::response::fillFileViews(input,products[0]->profile,classes,grids,weights,views);
 }
 photonjet::progressOffset+=photonjet::progressRowTotal;photonjet::progressPhase="selection_full_qa";
 std::vector<photonjet::qa::SelectionView> qaViews;for(size_t v=1;v<products.size();++v)if(!isJes(v))qaViews.push_back({products[v]->profile,&products[v]->qa,recipes[v],expectedGap(v)});
 photonjet::qa::fillFile(input,products[0]->profile,family,weights,classes,products[0]->qa,&recipe,"",{},qaViews,sim?nullptr:&mbdReference);
 photonjet::progressOffset+=photonjet::progressRowTotal;photonjet::progressPhase="selection_extensions";
 photonjet::fillExtensions(input.c_str(),(dirs[0]/"extensions.root").c_str(),(works[0]/"selection.json").c_str(),listed(sample.at("weights")).c_str(),h.at("require_stage3_fractions"),recipe.centralityOff);
 for(size_t v=0;v<products.size();++v){
  auto& p=*products[v];p.profile.simulation=sim;
  auto& measurement=isJes(v)?jesProducts.at(v)->histograms:p.measurement;
  auto& records=isJes(v)?jesProducts.at(v)->records:p.records;
  auto& timing=isJes(v)?jesProducts.at(v)->timing:p.timing;
  auto& timed=isJes(v)?jesProducts.at(v)->timingRecords:p.timed;
  {TFile f((dirs[v]/"spectra.root").c_str(),"CREATE");measurement.write();records.write();photonjet::writeSelection(p.profile);TObjString(p.profile.name.c_str()).Write("profile");
   if(!sim){f.mkdir("timecut")->cd();timing.write();timed.write();f.cd();if(system==photonjet::System::auau){const std::string trigger=plans[v].at("selection").at("systems").at("auau").at("data_trigger");f.mkdir(("triggers/"+trigger).c_str(),"",true);f.GetDirectory(("triggers/"+trigger).c_str())->cd();measurement.write();records.write();f.cd();TObjString(trigger.c_str()).Write("data_trigger");}}
  }
  if(isJes(v))fs::copy_file(dirs[0]/"qa.root",dirs[v]/"qa.root");else {TFile f((dirs[v]/"qa.root").c_str(),"CREATE");p.qa.write();photonjet::writeSelection(p.profile);TObjString(p.profile.name.c_str()).Write("profile");TObjString(photonjet::qa::familyName(family)).Write("family");
   TObjString(Json({{"schema","RJPhotonQAV1"},{"version",photonjet::qa::kQaVersion},{"et",{photonjet::qa::kQaEtMin,photonjet::qa::kQaEtMax}},{"stages",{"every candidate","preselected","tight","tight and isolated"}},{"histograms",p.qa.names().size()}}).dump().c_str()).Write("qa_contract");
   TObjString(Json({{"per_slice",weights.perSlice}}).dump().c_str()).Write("centrality_weights");TObjString(sourceText.c_str()).Write("sources");
   if(!sim){TFile in(input.c_str(),"READ");const char* key=in.Get("meta/trigger_selection")?"meta/trigger_selection":"meta/trigger_run_info";const auto info=photonjet::metaText(in,key);f.cd();TObjString(info.c_str()).Write("trigger_run_info");TObjString("1").Write("trigger_run_info_consistent");}
  }
  if(sim && (v || gridReceipt.is_null()))photonjet::response::writeProduct((dirs[v]/"response.root").c_str(),p.profile,grids,weights,recipes[v],p.response,sourceText);
  if(v){fs::copy_file(dirs[0]/"extensions.root",dirs[v]/"extensions.root");TFile f((dirs[v]/"extensions.root").c_str(),"UPDATE");photonjet::writeSelection(p.profile);TObjString(p.profile.name.c_str()).Write("profile",TObject::kOverwrite);}
 }
 Json denominator;
 if(sample.at("normalization").contains("ownership")){photonjet::progressOffset+=photonjet::progressRowTotal;photonjet::progressPhase="selection_stitching_denominator";denominator=stitchingDenominator(input,sample.at("normalization"));}
 require(stamp==std::make_pair(fs::file_size(input),fs::last_write_time(input)),"source changed during native fanout");
 for(size_t v=0;v<products.size();++v){
  const auto& p=plans[v];const auto code=digest(p.at("producer_code"));Json saved={{"input",input},{"sha256",source.at("expected_sha256")},{"root_version",gROOT->GetVersion()},{"plan_sha256",p.at("plan_sha256")},{"producer_code_sha256",code},{"range_begin",photonjet::entryBegin},{"range_end",photonjet::entryEnd},{"processed_entries",photonjet::progressRowTotal},{"sample",sample.at("key")},{"sample_class",sample.at("system").get<std::string>()+"/"+sample.at("family").get<std::string>()},{"reduction_settings_sha256",digest(Json{{"selection",p.at("selection")},{"histograms",p.at("histograms")}})},{"products",Json::object()},{"native_selection_fanout",{{"views",products.size()},{"physical_source_passes",sample.at("normalization").contains("ownership")?4:3},{"measurement_response_passes",1},{"full_qa_passes",1},{"cut_invariant_extensions_passes",1},{"extensions_rebound_from_parent",v!=0}}}};
  for(const std::string kind:sim?std::vector<std::string>{"spectra","qa","extensions","response"}:std::vector<std::string>{"spectra","qa","extensions"}){
   const auto hash=sha256File((dirs[v]/(kind+".root")).string());saved["products"][kind]=hash;
   writeText(dirs[v]/(kind+".json"),Json{{"plan_sha256",p.at("plan_sha256")},{"input_sha256",source.at("expected_sha256")},{"producer_code_sha256",code},{"range_begin",photonjet::entryBegin},{"range_end",photonjet::entryEnd},{"root_version",gROOT->GetVersion()},{"sha256",hash}}.dump()+"\n");
  }
  if(v==0 && !gridReceipt.is_null()){
   saved["native_response_grid"]=gridReceipt;
   for(const auto& product:gridReceipt.at("products"))if(product.at("name")!="nominal")saved["products"]["response__"+product.at("name").get<std::string>()]=product.at("sha256");
  }
  if(!denominator.is_null())saved["stitching_denominator"]=denominator;
  writeText(dirs[v]/"done.json",saved.dump(1)+"\n");
 }
}
void group(const std::vector<fs::path>& homes,const std::string& id,const fs::path& scratch){
 require(std::regex_match(id,std::regex("G[0-9]+")) && homes.size()>1 && homes.size()<=16,"invalid selection group or view count");require(fs::is_directory(scratch),"existing scratch directory required");
 const std::string shard=std::to_string(std::stoul(id.substr(1))/1000);std::vector<Json> descriptors,plans;std::set<fs::path> distinct;
 for(const auto& home:homes){require(distinct.insert(fs::canonical(home)).second,"duplicate view work");descriptors.push_back(groupPlan(home,id));plans.push_back(descriptors.back().at("plan"));checkReductionCode(plans.back());}
 const auto& parent=plans[0];require(!parent.contains("selection_fanout"),"first work must be the unchanged reference");
 require(!parent.at("histograms").value("candidate_records",false),"candidate archive fanout is not qualified");
 for(size_t v=1;v<plans.size();++v){const auto& p=plans[v];auto contract=p.at("selection_fanout");require(contract.at("parent_plan_sha256")==parent.at("plan_sha256") && contract.at("parent_selection_sha256")==parent.at("selection").at("selection_sha256"),"view reference differs");
  auto a=p,b=parent;a.erase("selection");b.erase("selection");b.erase("response_component_grid");if(contract.value("kind","")=="DATA_JES"){
   require(parent.at("samples").at(0).at("system")=="PP" && parent.at("samples").at(0).at("family")=="data","JES namespace on non pp DATA");
   const auto name=contract.at("prescription").get<std::string>();const double scale=name=="data_jes_up"?1./.907:name=="data_jes_down"?1./.923:0.;
   require(scale && a.at("histograms").at("data_jet_scale").at("pp")==scale && p.at("selection")==parent.at("selection"),"JES scale/selection differs");
   a["histograms"]["data_jet_scale"]["pp"]=b.at("histograms").at("data_jet_scale").at("pp");
  }
  if(parent.contains("response_component_grid")){
   require(a.at("histograms").at("response_variation")==viewRecipe(parent) && contract.at("reference_response_variation")==viewRecipe(parent),"cut view response differs from declared smnom component");
   a["histograms"]["response_variation"]=b.at("histograms").at("response_variation");
  }
  a.erase("selection_fanout");a.erase("plan_sha256");b.erase("plan_sha256");require(a==b,"group differs beyond the declared cuts");
  for(const auto& pin:contract.at("definition_sources"))require(sha256File(pin.at("path"))==pin.at("sha256").get<std::string>(),"view definition changed");
 }
 for(size_t v=0;v<homes.size();++v)require(!fs::exists(homes[v]/"chunks"/shard/id/"done.json")&&!fs::exists(homes[v]/"chunks"/shard/id/(parent.at("samples").at(0).at("key").get<std::string>()+".root")),"existing group output: reuse or reconcile before another traversal");
 const fs::path temporary=fs::canonical(scratch)/("selection-"+id+"-"+digest(Json{{"time",std::chrono::steady_clock::now().time_since_epoch().count()}}).substr(0,16));fs::create_directory(temporary);
 std::vector<fs::path> local;for(size_t v=0;v<homes.size();++v){local.push_back(temporary/std::to_string(v));savePlan(plans[v],local.back().string());}
 long long total=0,completed=0;for(const auto& j:parent.at("jobs"))total+=(j.at("range_end").get<long long>()-j.at("range_begin").get<long long>())*(parent.at("samples").at(0).at("normalization").contains("ownership")?4:3);
 std::map<std::string,double> timingReferences;
 if(parent.at("samples").at(0).at("family")=="data")for(const auto& source:parent.at("samples").at(0).at("sources"))total+=source.at("metadata").at("entries").get<long long>();
 const auto started=std::chrono::steady_clock::now();const double deadlineSeconds=parent.contains("response_component_grid")?parent.at("response_component_grid").at("limits").at("max_seconds").get<double>():28800.;
 const auto progress=[&](const std::string& phase,long long processed,const Json& error=nullptr){const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();require(phase=="FIRST_BAD" || seconds<=deadlineSeconds,"selection group runtime deadline exceeded");require(processed>=0 && processed<=total,"selection group progress exceeds declared work");const double rate=seconds>0?processed/seconds:0.;Json value={{"task","THE-510"},{"workstream","selection_fanout"},{"campaign",parent.at("plan_sha256")},{"row",id},{"phase",phase},{"processed",processed},{"total",total},{"remaining",total-processed},{"row_processed",processed},{"row_total",total},{"row_remaining",total-processed},{"deterministic_skip_processed",0},{"deterministic_skip_total",0},{"deterministic_skip_remaining",0},{"rate",rate},{"eta_seconds",rate>0?(total-processed)/rate:-1.},{"updated_at_unix",std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count()},{"latest_error",error},{"views",homes.size()},{"unit","physical_source_event_passes"}};writeText(temporary/"progress.json",value.dump()+"\n");};
 std::cout<<Json{{"schema","PhotonJetNativeSelectionGroupStartedV2"},{"group",id},{"progress_path",(temporary/"progress.json").string()},{"receipt_path",(temporary/"RECEIPT.json").string()},{"physical_source_event_passes",total}}.dump()<<std::endl;
 try{
  progress("starting",0);
  if(parent.at("samples").at(0).at("family")=="data")for(const auto& source:parent.at("samples").at(0).at("sources")){
    const std::string path=source.at("path");TFile file(path.c_str(),"READ");
    timingReferences[path]=photonjet::qa::mbdTimeReference(file,[&](Long64_t n){if(n%1024==0 || n==source.at("metadata").at("entries").get<Long64_t>())progress("source_timing_reference",completed+n);});
    completed+=source.at("metadata").at("entries").get<long long>();
  }
  photonjet::progressObserver=[&](const Json& child){progress(child.value("phase","reading"),completed+child.at("processed").get<long long>());};
  for(size_t i=0;i<parent.at("jobs").size();++i){row(plans,local,i,timingReferences.empty()?0.:timingReferences.at(parent.at("jobs").at(i).at("input").get<std::string>()));completed+=photonjet::progressTotal;}
  photonjet::progressObserver={};progress("merging",completed);
  for(size_t v=0;v<plans.size();++v){const auto& p=plans[v];const auto& sample=p.at("samples").at(0);merge(p,sample,local[v].string(),(local[v]/"packages").string());const auto original=local[v]/"packages"/(sample.at("key").get<std::string>()+".root");
   const auto output=homes[v]/"chunks"/shard/id;fs::create_directories(output);const auto target=output/(sample.at("key").get<std::string>()+".root");fs::create_directories(target.parent_path());const auto partial=target.string()+".selection-partial";fs::copy_file(original,partial);require(sha256File(partial)==sha256File(original.string()),"publication bytes differ");require(::link(partial.c_str(),target.c_str())==0,"publication collision requires reconciliation");fs::remove(partial);validatePackage(target.string(),p,sample);
   long long entries=0;std::set<std::string> processedSources;for(const auto& j:p.at("jobs")){entries+=j.at("range_end").get<long long>()-j.at("range_begin").get<long long>();processedSources.insert(j.at("input").get<std::string>());}
   const Json binding={{"plan_sha256",p.at("plan_sha256")},{"descriptor_sha256",descriptors[v].at("descriptor_sha256")},{"sha256",sha256File(target.string())}};writeText(output/"publication.json",binding.dump()+"\n");Json done=binding;done.update({{"status","COMPLETE"},{"processed_entries",entries},{"source_count",processedSources.size()},{"sample_class",sample.at("system").get<std::string>()+"/"+sample.at("family").get<std::string>()},{"seconds",nullptr},{"root_version",gROOT->GetVersion()},{"producer_code_sha256",digest(p.at("producer_code"))},{"reduction_settings_sha256",digest(Json{{"selection",p.at("selection")},{"histograms",p.at("histograms")}})},{"native_selection_fanout",{{"physical_work_owner",parent.at("plan_sha256")},{"views",homes.size()},{"science_acceptance",false}}}});writeText(output/"done.json",done.dump(1)+"\n");
  }
  progress("complete",total);writeText(temporary/"RECEIPT.json",Json{{"status","PASS_NATIVE_SELECTION_GROUP"},{"views",homes.size()},{"physical_source_event_passes",total},{"science_acceptance",false}}.dump(2)+"\n");
 }catch(const std::exception& e){photonjet::progressObserver={};try{progress("FIRST_BAD",completed,e.what());}catch(...){}throw;}
}
Json validateSaved(const std::vector<fs::path>& homes,const std::string& id){
 const std::string shard=std::to_string(std::stoul(id.substr(1))/1000);Json products=Json::array();std::set<std::string> namespaces;
 require(homes.size()>1 && std::regex_match(id,std::regex("G[0-9]+")),"invalid saved group/view roster");
 for(const auto& home:homes){const auto descriptor=groupPlan(home,id);const auto& p=descriptor.at("plan");checkReductionCode(p);const auto& sample=p.at("samples").at(0);
  const auto directory=home/"chunks"/shard/id;const auto target=directory/(sample.at("key").get<std::string>()+".root");
  const auto done=readJson((directory/"done.json").string());const auto publication=readJson((directory/"publication.json").string());const auto hash=sha256File(target.string());
  require(done.at("status")=="COMPLETE" && done.at("sha256")==hash && publication.at("sha256")==hash && done.at("descriptor_sha256")==descriptor.at("descriptor_sha256") && publication.at("descriptor_sha256")==descriptor.at("descriptor_sha256") && done.at("plan_sha256")==p.at("plan_sha256") && publication.at("plan_sha256")==p.at("plan_sha256"),"saved package binding differs");
  validatePackage(target.string(),p,sample);
  const std::string viewName=p.contains("selection_fanout")?p.at("selection_fanout").at("prescription").get<std::string>():"reference";require(namespaces.insert(viewName).second,"duplicate saved view namespace");
  products.push_back({{"namespace",viewName},{"package",fs::absolute(target).string()},{"sha256",hash},{"plan_sha256",p.at("plan_sha256")},{"selection_sha256",p.at("selection").at("selection_sha256")},{"components",p.value("response_component_grid",Json::object())},{"scope","COMPLETE_ORDINARY_PACKAGE"}});
 }
 return {{"schema","PhotonJetSavedSelectionGroupValidationV2"},{"status","PASS_STRUCTURAL_ONLY"},{"group",id},{"products",products},{"science_acceptance",false}};
}
} // selection_fanout
} // treetohists
