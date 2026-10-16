"""Typed scientific witnesses for existing PhotonID producers.

No training, reduction or new calculation occurs here. Effective parameters,
loaded model numbers and producer-written input witnesses enter identity;
filenames and digests remain provenance. Old populations stay unknown.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
import struct

CATEGORIES=('selection','weights','normalization','models','projection')


def unknown(reason):return {'state':'unknown','reason':reason}
def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)
def portable(value):
    if isinstance(value,float):
        if not math.isfinite(value):raise ValueError('nonfinite scientific definition')
        return 'f64:'+struct.pack('>d',value if value else 0.).hex()
    if isinstance(value,dict):return {k:portable(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [portable(v) for v in value]
    return value


def at(value,path):
    for name in path.split('/'):
        if not isinstance(value,dict) or name not in value:return unknown('upstream producer did not record '+path)
        value=value[name]
    return value


def source_definitions(metadata):
    chunks=metadata.get('consolidation',{}).get('chunks')
    rows=([chunk.get('original_meta',{}).get('scientific_definitions',unknown('source chunk lacks native definitions'))
           for chunk in chunks if chunk['entry_end']>chunk['entry_begin']] if chunks is not None else
          [metadata.get('scientific_definitions',unknown('source lacks native definitions'))])
    if not rows:rows=[unknown('no native source calculation observed')]
    return sorted({canonical(row):row for row in rows}.values(),key=canonical)


def definition(population,recipe,observable,*,weighted=False,model=None,selection=None,extra=None):
    receipt=population.receipt if hasattr(population,'receipt') else population
    rows=receipt.get('scientific_inputs',[unknown('prepared population has no native input witness')])
    selected=[]
    for row in rows:
        fields={key:at(row,key) for key in ('event_acceptance','candidate_geometry','truth_classification','vertex')}
        feature={}
        for name in recipe['features']:
            if name in ('cluster_Et','cluster_Eta'):value=at(row,'candidate_geometry/kinematics')
            elif name=='vertexz':value=at(row,'vertex')
            elif name=='centrality':value=at(row,'event_acceptance/centrality_measurement')
            else:value=at(row,'shower_feature_views/H70')
            feature[name]=value
        fields['features']=feature;selected.append(fields)
    selected=sorted({canonical(row):row for row in selected}.values(),key=canonical)
    terms=[]
    def add(category,key,value):terms.append({'category':category,'id':key,'value':portable(value),'windows':[]})
    add('selection','upstream_input_definitions',selected)
    add('selection','population',{'domain':recipe.get('training_domain',recipe['domain']),
        'feature_order':recipe['features'],'feature_precision':'binary32',
        'event_flags':{'record':0,'nodes_ready':1,'vertex_valid':1,'terminal':0},
        'labels':'Inputs.labels/v1: signal source-role and explicit nonprompt background',
        'background_flattening':receipt.get('background_flattening'),
        'background_flattening_parameters':recipe.get('training',{}).get('background_flattening'),
        'sample_roster':recipe.get('sources'),
        'minimum_bias':{key:recipe.get('training',{}).get(key) for key in ('minimum_bias_branch','minimum_bias_pass')}})
    if selection is not None:add('selection','diagnostic_subset',selection)
    add('weights','learner_weights',receipt.get('scientific_weights',unknown('prepared learner weights lack effective spline witness')) if weighted else {'rule':'unit'})
    add('normalization','diagnostic',{'kind':'raw native sufficient statistics; plot normalization is recorded by its consumer'})
    add('models','loaded_function',model if model is not None else {'state':'not_applicable','reason':'this diagnostic does not evaluate a learned function'})
    add('projection','observable',{'method':observable,'parameters':extra or {}})
    return {'schema':'EffectiveObservableDefinitionV1','observable':observable,
        'channel':{'system':recipe['system'],'population':'prepared_photonid'},'terms':terms,
        'categories':{key:{'state':'measured','value':'explicit diagnostic term roster'} for key in CATEGORIES}}


_ready=False

def native():
    global _ready
    import ROOT
    if not _ready:
        base=Path(__file__).resolve().parents[2]
        ROOT.gSystem.Load('libTMVA');ROOT.gSystem.Load('libTMVAUtils')
        code=('#include '+json.dumps(str(base/'FinalAnalysis/src/RootScientificTrace.h'))+'\n'+
              '#include '+json.dumps(str(base/'TreeToHists/src/NativeModelDefinition.h'))+'\n'+r'''
namespace photonid_scientific {
void attach(TObject* object,const char* definition) {
  photonjet_science::attachDefinition(object,photonjet_science::Json::parse(definition));
  photonjet_science::setTag(object,"physics_covariance",photonjet_science::covariance("not_computed",{},0,"native diagnostic sufficient statistics; no full cross-bin covariance computed"));
}
std::string trace(const char* definition) {return photonjet_science::Trace::source(photonjet_science::Json::parse(definition)).json().dump();}
void columns(TObject* object,const char* definitions) {
  const auto source=photonjet_science::Json::parse(definitions);auto out=photonjet_science::Json::object();
  for(const auto& item:source.items())out[item.key()]=photonjet_science::Trace::source(item.value()).json();
  photonjet_science::setTag(object,"physics_column_traces",out);
}
std::string model(const char* path,const char* key,const char* features) {
  return photonjet_science::readRBDTDefinition(path,key,photonjet_science::Json::parse(features)).dump();
}
}
''')
        if not ROOT.gInterpreter.Declare(code):raise ValueError('native scientific metadata bridge did not compile')
        _ready=True
    return ROOT


def attach(object,value):native().photonid_scientific.attach(object,canonical(value))
def attach_columns(object,values):native().photonid_scientific.columns(object,canonical(values))
def record(directory,name,value):
    ROOT=native();directory.cd()
    ROOT.TObjString(canonical({'trace':json.loads(str(ROOT.photonid_scientific.trace(canonical(value))))})).Write('physics_'+name)


def calibration_definition(population,observable,*,recipe=None,model=None,selection=None,extra=None):
    if recipe is None:recipe=population.receipt.get('scientific_recipe') if hasattr(population,'receipt') else None
    if recipe is None:
        recipe={'system':'unknown','features':[],'domain':unknown('calibration population recipe not recorded')}
    return definition(population,recipe,observable,weighted=True,
        model=model if model is not None else unknown('score producer did not supply a numerical model witness'),
        selection=selection,extra=extra)


def numpy_function(artifact):
    features=artifact['features'];pre=artifact['preprocessing'];layers=artifact['layers']
    if not features or len(set(features))!=len(features) or len(pre['mean'])!=len(features) or len(pre['scale'])!=len(features):
        raise ValueError('numpy model feature/transform dimensions differ')
    if any(not math.isfinite(v) or v<=0 for v in pre['scale']):raise ValueError('numpy model scale is invalid')
    width=len(features)
    for layer in layers:
        if not layer['weights'] or len(layer['bias'])!=len(layer['weights']) or any(len(row)!=width for row in layer['weights']):
            raise ValueError('numpy model layer dimensions differ')
        width=len(layer['weights'])
    if width!=1 or artifact.get('hidden_activation')!='tanh' or artifact.get('output')!='sigmoid':
        raise ValueError('numpy model output or activation unsupported')
    return portable({'schema':'NumpyScientificFunctionV1','feature_order':features,'input_precision':'binary64',
        'preprocessing':pre,'layers':layers,'hidden_activation':'tanh',
        'output':'sign-stable sigmoid: x>=0 -> 1/(1+exp(-x)); else exp(x)/(1+exp(x))'})


def model_function(folder,recipe,cells):
    from . import Models
    result=[];folder=Path(folder).resolve()
    for cell in cells:
        path=(folder/cell['artifact']).resolve()
        Models.require(path.is_relative_to(folder),'scientific model escapes release')
        before=Models.file_hash(path);artifact=Models.read_json(path)
        Models.require(artifact['features']==recipe['features'],'scientific model feature order differs')
        if artifact['schema']=='PhotonIDNumpyModelV1':function=numpy_function(artifact)
        elif artifact['schema']=='PhotonIDXGBoostModelV1' and 'native_evaluator' in artifact:
            spec=artifact['native_evaluator'];product=(path.parent/spec['file']).resolve()
            Models.require(product.is_relative_to(folder),'scientific native model escapes release')
            digest=Models.file_hash(product)
            function=json.loads(str(native().photonid_scientific.model(str(product),spec['key'],canonical(recipe['features']))))
            Models.require(Models.file_hash(product)==digest,'native model changed during scientific readback')
        else:function=unknown('saved model has no supported numeric evaluator representation')
        Models.require(Models.file_hash(path)==before,'model artifact changed during scientific readback')
        result.append({'bounds':cell.get('bounds',{}),'function':function})
    return {'routing':'half-open cell bounds in declared order','domain':recipe['domain'],'cells':result}
