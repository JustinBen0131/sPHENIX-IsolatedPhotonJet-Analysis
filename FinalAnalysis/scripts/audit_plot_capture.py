#!/usr/bin/env python3
"""Verify an add-only legacy capture overlay preserves all original C++ tokens.

Only marked, allowlisted metadata calls can be added. Any changed calculation,
fit option, parameter, control flow, literal or source include refuses PASS.
This is deliberately narrower than a general C++ semantic equivalence proof.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
from plot_physics_inventory import sha,require

BEGIN='/* PHYSICS_CAPTURE_BEGIN */'
END='/* PHYSICS_CAPTURE_END */'
CALLS={'bind','bindNative','annotation','attachTrace','attachDefinition','setTag','decoration',
       'Trace::source','Trace::fromJson','Trace::combine','traceOf','projectionOf','selectedTrace',
       'nativeValues','histogram','bins','selectedHistogram','covarianceOf','scalar','covariance','real','known','absent','unknown'}
LITERALS=re.compile(r'R"([^ ()\\\t\r\n]*)\([\s\S]*?\)\1"|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*[\s\S]*?\*/')


def tokens(text):
    parts=[];at=0
    for match in LITERALS.finditer(text):
        parts.extend(re.findall(r'\w+|[^\w\s]',text[at:match.start()]))
        value=match.group()
        if not value.startswith(('//','/*')):parts.append(value)
        at=match.end()
    parts.extend(re.findall(r'\w+|[^\w\s]',text[at:]));return parts


def capture_only(block):
    if re.fullmatch(r'\s*#include "(?:\.\./)*(?:src/)?RootScientificTrace\.h"\s*',block):return
    clean=LITERALS.sub(' ',block)
    require(not re.search(r'(?<![=!<>])=(?!=)|\+\+|--|\b(?:new|delete|return|throw|while|for|if|goto|typedef|using)\b',clean),
            'capture block changes control flow, allocates or assigns')
    require(not re.search(r'\[[^\]]*\]\s*(?:\([^)]*\))?\s*\{',clean),'capture block contains a lambda')
    require(not re.search(r'[A-Za-z_][\w:]*\s*\{',clean),'capture block constructs an unreviewed type')
    statements=[];start=0;paren=brace=0
    for index,char in enumerate(clean):
        paren+=(char=='(')-(char==')');brace+=(char=='{')-(char=='}')
        require(paren>=0 and brace>=0,'unbalanced capture expression')
        if char==';' and paren==0 and brace==0:statements.append(clean[start:index].strip());start=index+1
    require(not clean[start:].strip() and paren==brace==0,'capture block must contain complete metadata statements')
    require(statements and all(re.match(r'^photonjet_science::[\w:]+\s*\(',line) and line.endswith(')') for line in statements),
            'capture block contains a declaration or non-capture statement')
    calls=re.findall(r'([A-Za-z_][\w:]*)\s*\(',clean)
    require(calls,'capture block has no capture call')
    for call in calls:
        allowed=(call.startswith('photonjet_science::') and call[len('photonjet_science::'):] in CALLS)
        allowed=allowed or call in {'at','get','value','then','json','GetXmin','GetXmax','GetNbinsX','GetNbinsY','GetParameter','GetParError','GetNpar','GetNDF','GetChisquare','GetBinContent','GetBinError','GetXaxis','GetYaxis','GetBinLowEdge','GetBinUpEdge','__builtin_FUNCTION'}
        require(allowed,'capture block calls an unreviewed operation: '+call)
    require(clean.lstrip().startswith('photonjet_science::'),'capture statement must be an explicit metadata API call')
    # Calls may contain literal JSON and getter arguments. Member assignment,
    # comma expressions and arbitrary helper evaluation are rejected above.


def strip_capture(text):
    require(text.count(BEGIN)==text.count(END),'unbalanced capture markers')
    result=[];at=0;count=0
    while BEGIN in text[at:]:
        begin=text.index(BEGIN,at);end=text.find(END,begin+len(BEGIN))
        require(end>=0 and BEGIN not in text[begin+len(BEGIN):end],'nested capture markers')
        capture_only(text[begin+len(BEGIN):end]);result.append(text[at:begin]);at=end+len(END);count+=1
    result.append(text[at:]);return ''.join(result),count


def audit(original_files,capture_files,original_pins,capture_pins):
    require(set(original_files)==set(original_pins),'old source roster differs from its exact receipt')
    require(set(original_files)<=set(capture_files),'capture lost an original source')
    evidence=[]
    for name,path in original_files.items():
        require(sha(path)==original_pins[name],'retained original source changed: '+name)
        require(sha(capture_files[name])==capture_pins[name],'capture source changed: '+name)
        original=Path(path).read_text();captured=Path(capture_files[name]).read_text()
        require(BEGIN not in original and END not in original,'original contains capture markers already')
        stripped,count=strip_capture(captured)
        require(tokens(original)==tokens(stripped),'a non-capture scientific/source token changed: '+name)
        evidence.append({'name':name,'original_sha256':sha(path),'capture_sha256':sha(capture_files[name]),'capture_blocks':count,
                         'preserved_token_sha256':hashlib.sha256(json.dumps(tokens(original)).encode()).hexdigest()})
    return {'schema':'PlotCaptureInstrumentationAuditV1','scientific_expressions_unchanged':True,
            'method':'add_only_allowlisted_capture_blocks_v1','original_renderer':original_pins,
            'capture_renderer':capture_pins,'source_token_evidence':evidence}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();m=json.loads(a.manifest.read_text());result=audit(m['original_files'],m['capture_files'],m['original_renderer'],m['capture_renderer'])
    with a.output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
