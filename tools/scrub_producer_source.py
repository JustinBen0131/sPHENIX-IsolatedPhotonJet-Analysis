#!/usr/bin/env python3
"""Deterministic clean-up of the producer sources copied into TreeProduction/.

The Fun4All producer (TreeProduction/src, src_AuAu, macros) is copied verbatim
from the production source at each production freeze and then passed through
this script.  The script is idempotent: running it on already-scrubbed sources
changes nothing, and ``--check`` verifies the result.

What it changes
  1. Internal study numbers (``THE-105``, ``THE-134``, ...) are removed from
     comments.  They are bookkeeping labels with no physics meaning.
  2. Absolute paths of one user's installation are replaced by the standard
     sPHENIX environment lookups or by an environment variable:
       * ``R__LOAD_LIBRARY(<abs>/libcalo_reco.so)`` -> ``R__LOAD_LIBRARY(libcalo_reco.so)``
       * ``#include "<abs>/caloreco/PhotonClusterBuilder.h"`` -> ``<caloreco/PhotonClusterBuilder.h>``
       * the scaled-trigger QA run list is taken from ``RJ_SCALED_TRIGGER_RUNLIST``
         (already honoured by the macro) instead of a hard-coded scratch path;
       * simulation reweighting files are configured only in analysis_config.yaml.

What it deliberately keeps
  * C++ identifiers, environment-variable names, ROOT histogram/tree names and
    log tags that still carry study numbers (``m_the44...``, ``RJ_THE134_*``,
    ``THE44PythiaAutopsyTree``, ``h_the119_*``).  Renaming them would change the
    runtime interface used by the production wrappers or the output files, and
    the producer cannot be compiled on the machine that runs this script.
  * ``THE106Observation``: the name of an optional instrumentation header in a
    patched PHOOL release; the local ``THE106ObservationDisabled.h`` mirrors its
    API and must keep the name.

Usage
  python tools/scrub_producer_source.py            # apply
  python tools/scrub_producer_source.py --check    # verify, exit 1 on findings
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
PRODUCER = ROOT / "TreeProduction"
SOURCE_SUFFIXES = {".C", ".cc", ".h"}

TASK_TOKEN = re.compile(r"\bTHE-?\d{2,3}(?:/\d{2,3})?\b")
# Paths of the producing author's own installation and scratch areas.  Input files
# provided by collaborators (calibration masks, reweighting histograms) are
# external inputs listed in TreeProduction/README.md and are allowed.
SITE_PATH = re.compile(r"/sphenix/(?:u|user)/(?:patsfan753|jbennett)/|/sphenix/tg/tg01/bulk/jbennett|/gpfs/")

# Phrases whose meaning would be garbled by plain token removal.
PHRASES = {
    "// THE-117 certifies that the authoritative direct path builds and runs without":
        "// The authoritative direct path is certified to build and run without",
}

# Exact, one-time source rewrites for site-specific paths.  Each entry is
# (relative file, old text, new text).  Missing "old" text is not an error so
# the script stays idempotent.
SITE_REWRITES = [
    (
        "macros/Calo_Calib.C",
        "R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)",
        "R__LOAD_LIBRARY(libcalo_reco.so)",
    ),
    (
        "macros/Fun4All_recoilJets_unified_impl.C",
        '#include "/sphenix/u/patsfan753/thesisAnalysis/install/include/caloreco/PhotonClusterBuilder.h"',
        "#include <caloreco/PhotonClusterBuilder.h>",
    ),
    (
        "macros/Fun4All_recoilJets_unified_impl.C",
        '  std::string m_runListPath =\n'
        '    "/sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau/"\n'
        '    "scaledEffRuns_MBD_NS_geq_2_vtx_lt_150__Pho10_12.list";',
        '  // Optional scaled-trigger QA run list; set RJ_SCALED_TRIGGER_RUNLIST to use it.\n'
        '  std::string m_runListPath;',
    ),
    (
        "macros/Fun4All_recoilJets_unified_impl.C",
        '        std::string vertex_reweight_file_pp = "/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/0mrad/reweight.root";',
        '        std::string vertex_reweight_file_pp;    // set in analysis_config.yaml',
    ),
    (
        "macros/Fun4All_recoilJets_unified_impl.C",
        '        std::string vertex_reweight_file_auau = "/sphenix/u/bseidlitz/work/pj_auau/reweightingDer/output/vtxz_reweighting.root";',
        '        std::string vertex_reweight_file_auau;  // set in analysis_config.yaml',
    ),
    (
        "macros/Fun4All_recoilJets_unified_impl.C",
        '        std::string centrality_reweight_file = "/sphenix/u/bseidlitz/work/pj_auau/reweightingDer/output/centrality_reweighting.root";',
        '        std::string centrality_reweight_file;   // set in analysis_config.yaml',
    ),
    (
        "src_AuAu/RecoilJets_AuAu.cc",
        '                    const std::string scaledTrigRunListPath =\n'
        '                    "/sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau/"\n'
        '                    "scaledEffRuns_MBD_NS_geq_2_vtx_lt_150__Pho10_12.list";',
        '                    // Optional scaled-trigger QA run list, shared with the steering macro.\n'
        '                    const char* scaledTrigRunListEnv = std::getenv("RJ_SCALED_TRIGGER_RUNLIST");\n'
        '                    const std::string scaledTrigRunListPath =\n'
        '                    scaledTrigRunListEnv ? std::string(scaledTrigRunListEnv) : std::string();',
    ),
]


def source_files() -> list[Path]:
    return sorted(p for p in PRODUCER.rglob("*") if p.is_file() and p.suffix in SOURCE_SUFFIXES)


def scrub_comment_text(text: str) -> str:
    """Remove task tokens from one comment fragment and tidy the spacing."""

    cleaned = TASK_TOKEN.sub("", text)
    cleaned = re.sub(r"[ \t]{2,}", " ", cleaned)
    cleaned = re.sub(r"\(\s*\)", "", cleaned)
    cleaned = re.sub(r"//\s+", "// ", cleaned, count=1) if cleaned.lstrip().startswith("//") else cleaned
    return cleaned.rstrip() if text.endswith(("\n", "")) else cleaned


def scrub_comments(text: str) -> str:
    """Apply phrase overrides, then remove task tokens inside comments only."""

    for old, new in PHRASES.items():
        text = text.replace(old, new)
    out_lines = []
    in_block = False
    for line in text.split("\n"):
        stripped = line.lstrip()
        if in_block:
            body = line
            if "*/" in line:
                in_block = False
        elif stripped.startswith("/*"):
            in_block = "*/" not in line
            body = line
        elif "//" in line:
            code, comment = line.split("//", 1)
            if TASK_TOKEN.search(comment):
                comment = TASK_TOKEN.sub("", comment)
                comment = re.sub(r"[ \t]{2,}", " ", comment)
                comment = re.sub(r"\(\s*\)", "", comment)
                # Capitalise a comment that now starts with a lowercase word.
                m = re.match(r"(\s*)([a-z])", comment)
                if m and code.strip() == "":
                    comment = m.group(1) + m.group(2).upper() + comment[m.end():]
            out_lines.append(code + "//" + comment.rstrip() if TASK_TOKEN.search(line) else line)
            continue
        else:
            out_lines.append(line)
            continue
        # block-comment line
        if TASK_TOKEN.search(body):
            body = TASK_TOKEN.sub("", body)
            body = re.sub(r"[ \t]{2,}", " ", body).rstrip()
        out_lines.append(body)
    return "\n".join(out_lines)


def apply() -> int:
    changed = 0
    for rel, old, new in SITE_REWRITES:
        path = PRODUCER / rel
        text = path.read_text(encoding="utf-8", errors="surrogateescape")
        if old in text:
            path.write_text(text.replace(old, new, 1), encoding="utf-8", errors="surrogateescape")
            changed += 1
            print(f"[scrub] site path: {rel}")
    for path in source_files():
        text = path.read_text(encoding="utf-8", errors="surrogateescape")
        scrubbed = scrub_comments(text)
        if scrubbed != text:
            path.write_text(scrubbed, encoding="utf-8", errors="surrogateescape")
            changed += 1
            print(f"[scrub] comments: {path.relative_to(ROOT)}")
    print(f"[scrub] {changed} file edit(s)")
    return 0


def check() -> int:
    findings: list[str] = []
    kept: dict[str, int] = {}
    for path in source_files():
        rel = path.relative_to(ROOT)
        for number, line in enumerate(path.read_text(encoding="utf-8", errors="surrogateescape").split("\n"), 1):
            stripped = line.lstrip()
            comment = None
            if stripped.startswith(("//", "*", "/*")):
                comment = line
            elif "//" in line:
                comment = line.split("//", 1)[1]
            if comment is not None and TASK_TOKEN.search(comment):
                findings.append(f"{rel}:{number}: task token in comment: {line.strip()[:100]}")
            if SITE_PATH.search(line) and not stripped.startswith(("//", "*")):
                findings.append(f"{rel}:{number}: site-specific absolute path: {line.strip()[:100]}")
            for token in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", line):
                if re.search(r"(?i)the\d{2,3}", token):
                    kept[token] = kept.get(token, 0) + 1
    yaml_path = PRODUCER / "macros" / "analysis_config.yaml"
    for number, line in enumerate(yaml_path.read_text(encoding="utf-8").split("\n"), 1):
        if SITE_PATH.search(line) and not line.lstrip().startswith("#") and "_file" not in line.split(":", 1)[0]:
            findings.append(f"{yaml_path.relative_to(ROOT)}:{number}: site path outside an external-input key: {line.strip()[:100]}")
    print(f"[check] identifiers with study numbers kept on purpose: {len(kept)} distinct")
    for line in findings:
        print("[check] " + line)
    print(f"[check] {'PASS' if not findings else 'FAIL'}: {len(findings)} finding(s)")
    return 0 if not findings else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="verify only; exit 1 on findings")
    args = parser.parse_args(argv)
    return check() if args.check else apply()


if __name__ == "__main__":
    sys.exit(main())
