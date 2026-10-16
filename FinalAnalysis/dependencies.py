#!/usr/bin/env python3
"""Local C++ dependencies for final_analysis.sh's content-addressed build.

Ask the selected compiler for the source's transitive includes with the same
flags as compilation. Hash only files inside the sibling-package source tree;
ROOT and the compiler retain their version identities in the launcher. No
configuration, model or data payload is read unless C++ actually includes it.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys


def local_pins(compiler, flags, package, sources):
    package = Path(package).resolve()
    root = package.parent
    pins = {}
    for source in sources:
        result = subprocess.run(
            [compiler, *shlex.split(flags), "-I" + str(package / "src"),
             "-I" + str(package / "plots"), "-MM", "-MT", "dependencies",
             "-x", "c++", source], check=True, text=True, stdout=subprocess.PIPE)
        # GCC/Clang emit Make dependencies, not shell words: quotes in a path
        # are literal, spaces/# are backslash escaped and dollar signs doubled.
        rule = re.sub(r"\\\r?\n", "", result.stdout)
        target, separator, prerequisites = rule.partition(":")
        if target != "dependencies" or not separator:
            raise ValueError("unrecognized compiler dependency output")
        for token in re.findall(r"(?:\\.|[^\s])+", prerequisites):
            name = re.sub(r"\\([\\ #:\t])", r"\1", token).replace("$$", "$")
            path = Path(os.path.abspath(name))
            try:
                relative = path.relative_to(root)
            except ValueError:
                continue
            pins[relative.as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pins


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.exit("usage: dependencies.py COMPILER FLAGS PACKAGE SOURCE...")
    try:
        print(json.dumps(local_pins(*sys.argv[1:4], sys.argv[4:]), sort_keys=True, separators=(",", ":")))
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        sys.exit("FinalAnalysis dependencies: " + str(error))
