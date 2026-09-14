#!/usr/bin/env python3
"""Add replay capture to an exact, already prepared collaborator source.

Writes a local overlay only. Explicit node bindings are required for every lane;
this tool cannot decide which embedded background nodes are scientifically valid.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from centrality_replay import SCALARS, ARRAYS, NODES, ALGORITHM


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"nonunique source anchor: {old[:80]}")
    return text.replace(old, new, 1)


def prepare(binding_path, output):
    binding = json.loads(binding_path.read_text())
    if binding.get("schema") != "CentralityCaptureSourceBindingV1":
        raise ValueError("unsupported source binding")
    nodes = binding["nodes"]
    for k in NODES:
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", nodes.get(k, "")):
            raise ValueError(f"invalid explicit node binding: {k}")
    if binding["sample_kind"] not in ("data", "embedded"):
        raise ValueError("capture overlay is for AuAu DATA or embedding")
    if binding["sample_kind"] == "embedded" and not binding.get("background_nodes_verified"):
        raise ValueError("embedded background node selection requires evidence")
    root = Path(binding["source_root"])
    header = "src/RJReplayFoundationV1.h"
    source = binding["emission_source"]
    if source not in ("src_AuAu/RecoilJets_AuAu.cc", "src_AuAu/detail/ReplayEventEmission.inc"):
        raise ValueError("unsupported emission source")
    contents = {}
    for name in (header, source):
        if digest(root / name) != binding["source_sha256"][name]:
            raise ValueError(f"source drift: {name}")
        contents[name] = (root / name).read_text()
    fields = "\n".join(f"  {'std::int32_t' if kind=='int32' else 'double'} {name}=" +
                       ("-1;" if kind == "int32" else "std::numeric_limits<double>::quiet_NaN();")
                       for name, kind in SCALARS.items())
    fields += "\n" + "\n".join(f"  std::vector<{'int' if kind.endswith('int32') else 'double'}> {name};"
                               for name, kind in ARRAYS.items())
    anchor = "  double mbd_total_charge=std::numeric_limits<double>::quiet_NaN();"
    contents[header] = once(contents[header], anchor, anchor + "\n" + fields)
    anchor = 'scalar(m_tEvent,"mbd_total_charge",&m_event.mbd_total_charge,"D");'
    bookings = "".join(f' scalar(m_tEvent,"{n}",&m_event.{n},"{"I" if t=="int32" else "D"}");'
                       for n, t in SCALARS.items())
    bookings += "".join(f' vec(m_tEvent,"{n}",&m_event.{n});' for n in ARRAYS)
    contents[header] = once(contents[header], anchor, anchor + bookings)
    # Capture at the beginning of event serialization, including terminal and
    # truth-only rows. The original event occurrence and run identities survive.
    anchor = "RJReplayRuntimeV1::EventBundle bundle;"
    call = f'''
    RJCentralityReplayV1::capture(bundle.event,
      findNode::getClass<MbdPmtContainer>(topNode,"{nodes['pmt']}"),
      findNode::getClass<MbdOut>(topNode,"{nodes['mbd_out']}"),
      findNode::getClass<MinimumBiasInfo>(topNode,"{nodes['minimum_bias']}"));
    // The same empty CentralityInfov2::Reset() makes the node stale for a
    // non-minimum-bias event.  Record the native witness only where
    // CentralityReco actually wrote it for this event; embedded input
    // samples carry their own per-event centrality product instead.
    const bool nativeCentralityTrusted=
      m_isSimEmbedded||bundle.event.centrality_mb_decision==1;
    auto* replayCentrality=nativeCentralityTrusted
      ? findNode::getClass<CentralityInfo>(topNode,"{nodes['centrality']}")
      : nullptr;
    bundle.event.centrality_native_valid=0;
    if(replayCentrality && replayCentrality->has_centile(CentralityInfo::PROP::mbd_NS)
       && replayCentrality->has_centrality_bin(CentralityInfo::PROP::mbd_NS))
    {{
      bundle.event.centrality_native_centile=replayCentrality->get_centile(CentralityInfo::PROP::mbd_NS);
      bundle.event.centrality_native_bin=replayCentrality->get_centrality_bin(CentralityInfo::PROP::mbd_NS);
      const auto c=bundle.event.centrality_native_centile;
      const auto b=bundle.event.centrality_native_bin;
      bundle.event.centrality_native_valid=(nativeCentralityTrusted && std::isfinite(c)
          && c>0 && c<=1 && b>=1 && b<=100) ? 1 : 0;
    }}
'''
    inc = "../../src/" if source.endswith(".inc") else "../src/"
    includes = (f'#include "{inc}RJCentralityReplayV1.h"\n#include <mbd/MbdPmtContainer.h>\n'
                '#include <mbd/MbdPmtHit.h>\n#include <mbd/MbdOut.h>\n'
                '#include <calotrigger/MinimumBiasInfo.h>\n#include <centrality/CentralityInfo.h>\n')
    contents[source] = includes + once(contents[source], anchor, anchor + call)
    contents["src/RJCentralityReplayV1.h"] = Path(__file__).with_name("RJCentralityReplayV1.h").read_text()
    output.mkdir(parents=True, exist_ok=False)
    for name, text in contents.items():
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)
    receipt = {"schema": "CentralityCaptureSourceDeltaV1", "status": "PREPARED_SOURCE_ONLY",
               "binding": str(binding_path.resolve()), "binding_sha256": digest(binding_path),
               "nodes": nodes, "algorithm": ALGORITHM,
               "outputs": {name: digest(output / name) for name in contents},
               "compiled_in_sphenix": False, "real_input_parity": False, "production_allowed": False}
    (output / "SOURCE_DELTA.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binding", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(prepare(a.binding, a.output), indent=2))
