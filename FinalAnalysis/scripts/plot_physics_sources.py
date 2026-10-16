#!/usr/bin/env python3
"""Inventory qualification sources and maintained native save entrypoints.

Source hashes qualify the instrumentation release. They are never scientific
values, and this static roster is never evidence that a family fixture ran.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re

PACKAGES=('FinalAnalysis','TreeToHists','TreeProduction','PhotonID')
SOURCE_SUFFIXES={'.C','.cc','.cpp','.cxx','.h','.hh','.hpp','.py','.json','.yaml','.yml','.sh','.txt','.cmake','.am','.ac'}
EXCLUDE={'build','outputs','output','runs','run','cache','__pycache__','.git','weights','fig','.physics'}
EXCLUDE_NAMES={'plot_physics_coverage.json','plot_physics_qualification.json','physics_family_roster.json',
               'BuildInfo.h','ReplayBuildInfo.h','BuildIdentity.h','EventTreeSchema.h'}
WRAPPERS={'save','savePhysicsQA','saveModelQA','saveReleasedWorkingPoints','saveTrainingDiagnostic'}


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_roots():
    """The exact same closed scan specification used by the canonical owner."""
    return [{'path':package,'suffixes':sorted(SOURCE_SUFFIXES),'names':['CMakeLists.txt','Makefile'],
             'exclude_dirs':sorted(EXCLUDE),'exclude_names':sorted(EXCLUDE_NAMES)} for package in PACKAGES]


def source_inventory(repository):
    root=Path(repository).resolve();result={}
    for spec in source_roots():
        directory=root/spec['path']
        for path in sorted(directory.rglob('*')):
            relative=path.relative_to(directory)
            if (not path.is_file() or set(relative.parts[:-1]) & set(spec['exclude_dirs']) or
                    path.name in spec['exclude_names'] or
                    (path.suffix not in spec['suffixes'] and path.name not in spec['names'])):continue
            resolved=path.resolve()
            if not resolved.is_relative_to(root):raise ValueError('source inventory escapes repository: '+str(path))
            result[path.relative_to(root).as_posix()]=digest(path)
    return result


def without_literals(text):
    # Preserve newlines and character positions for a bounded lexical call-site
    # inventory. This is not an AST equivalence checker or a science parser.
    pattern=r'R"([^ ()\\\t\r\n]*)\([\s\S]*?\)\1"|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*[\s\S]*?\*/'
    return re.sub(pattern,lambda m:''.join('\n' if c=='\n' else ' ' for c in m.group()),text)


def families(analysis):
    result={}
    for filename,prefix in [('IANFigures.C','ian.'),('QAFigures.C','qa.')]:
        path=Path(analysis)/'plots'/filename;text=path.read_text();clean=without_literals(text)
        header=re.compile(r'^[A-Za-z_][A-Za-z_0-9:<>, *&]*?\s+([A-Za-z_]\w*)\s*\([^;{}]*\)\s*\{',re.M)
        for match in header.finditer(clean):
            name=match.group(1)
            if name in WRAPPERS:continue
            at=match.end();depth=1
            while depth and at<len(clean):
                depth+=(clean[at]=='{')-(clean[at]=='}');at+=1
            body=clean[match.end():at]
            saves=sorted({call for call in WRAPPERS if re.search(r'\b'+call+r'\s*\(',body)})
            if not saves:continue
            key=prefix+name
            row={'source':'plots/'+filename,'entrypoint':name,'line':text.count('\n',0,match.start())+1,
                 'save_routes':saves,'qualification':'REQUIRES_EXECUTED_NATIVE_FAMILY_FIXTURE'}
            if key in result:
                row['overloads']=result[key].get('overloads',1)+1
            result[key]=row
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--repository',type=Path,required=True)
    args=parser.parse_args();root=args.repository.resolve()
    print(json.dumps({'schema':'PlotPhysicsSourceRosterV1','source_root':str(root),'source_roots':source_roots(),
                      'source_inventory':source_inventory(root),'families':families(root/'FinalAnalysis')},indent=2))

if __name__=='__main__':main()
