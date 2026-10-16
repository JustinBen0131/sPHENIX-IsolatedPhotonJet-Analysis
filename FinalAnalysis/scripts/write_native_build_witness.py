#!/usr/bin/env python3
"""Embed the complete producer source inventory BEFORE compiling a binary."""
import argparse
import json
import os
from pathlib import Path
import tempfile
import hashlib
from plot_physics_sources import source_inventory, source_roots


def write(source_root,output=None,output_dir=None):
    root=Path(source_root).resolve()
    sources=source_inventory(root)
    if not sources:raise ValueError('empty native producer source inventory')
    witness={'schema':'NativeProducerSourceWitnessV1','state':'recorded',
             'source_root':str(root),'source_roots':source_roots(),'source_inventory':sources}
    value=json.dumps(witness,sort_keys=True,separators=(',',':'),ensure_ascii=True)
    header='#pragma once\n#define PHOTONJET_NATIVE_SOURCE_WITNESS '+json.dumps(value,ensure_ascii=True)+'\n'
    if output_dir is not None:
        directory=Path(output_dir).resolve();directory.mkdir(parents=True,exist_ok=True)
        output=directory/(hashlib.sha256(header.encode()).hexdigest()+'.h')
        try:
            with output.open('x') as stream:stream.write(header);stream.flush();os.fsync(stream.fileno())
        except FileExistsError:
            if output.read_text()!=header:raise ValueError('immutable build witness collision')
        if source_inventory(root)!=sources:raise ValueError('producer source changed while creating build witness')
        return output
    output=Path(output).resolve()
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists() and output.read_text()==header:return witness
    with tempfile.NamedTemporaryFile('w',dir=output.parent,delete=False) as stream:
        stream.write(header);stream.flush();os.fsync(stream.fileno());name=stream.name
    if source_inventory(root)!=sources:raise ValueError('producer source changed while creating build witness')
    os.replace(name,output);return witness


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-root',type=Path,required=True)
    group=p.add_mutually_exclusive_group(required=True);group.add_argument('--output',type=Path);group.add_argument('--output-dir',type=Path)
    a=p.parse_args()
    value=write(a.source_root,a.output,a.output_dir)
    if a.output_dir:print(str(value))
    else:print(json.dumps({'sources':len(value['source_inventory']),'output':str(a.output.resolve())}))
if __name__=='__main__':main()
