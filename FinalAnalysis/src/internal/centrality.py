"""Au+Au simulation centrality weights: w(c_i) = p_DATA(c_i) / p_SIM(c_i) for the 5% slices c_i of 0-80%.

    python -m FinalAnalysis.src.internal.centrality --config config/final_analysis.yaml --data @auau_qa_data_outputs.txt \
        --simulation auau_photon12 @sums12.txt --simulation auau_photon20 @sums20.txt --out weights.json \
        [--data-trigger "Photon 12 GeV + MBD NS >= 2, vtx < 150 cm"]

p_X(c_i) is the fraction of the centrality distribution of X inside 0-80% that falls in slice i.
DATA: the events that fired the named scaled trigger, minimum-bias classified and inside the analysis vertex window, from
the FillPhotonQA data outputs (cent_bits_data, one run per output). The trigger is found by name in each run's own
trigger table, so bit numbers may differ between run eras, and a run where it is absent or disabled (initial prescale
< 0) contributes nothing and is listed. The default is scaled Photon 12 GeV +
MBD NS >= 2, vtx < 150 cm; the Photon 10 alternative uses its own DATA target and weights.
Simulation: the histogram of every simulation file (hist/SIM/h_centrality), summed per sample and then over the samples
with each sample's per-event normalisation (package.sample_weight), as the samples are combined everywhere else.
Slices are [5i, 5i + 5) on the histograms' own edges, the same floor(c / 5) FillPhotonJetResponse.C applies the weights
with. The output lists the 16 weights (the FillPhotonJetResponse.C argument) with the counts and the sha256 of every input.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import uproot
import yaml

from . import inputs as package
from .inputs import AUAU_DATA_TRIGGERS

DATA_TRIGGER = AUAU_DATA_TRIGGERS['photon12']
DATA_BITS = "cent_bits_data"                                 # FillPhotonQA.C: centrality x scaled bit
SIM_HISTOGRAM = "hist/SIM/h_centrality"
SLICE_EDGES = np.arange(0.0, 81.0, 5.0)


def summed(paths, key) -> tuple[np.ndarray, np.ndarray]:
    """(values, edges) of one histogram summed over files with identical edges."""
    total, edges = None, None
    for path in paths:
        with uproot.open(path) as f:
            if key not in f:
                raise ValueError(f"{path} has no {key}: its production lane does not carry the centrality histogram")
            h = f[key]
            values, these = h.values(), h.axes[0].edges()
        if edges is None:
            total, edges = values.astype(float).copy(), these
        elif not np.array_equal(these, edges):
            raise ValueError(f"{path}: {key} has other bin edges")
        else:
            total += values
    if total is None:
        raise ValueError(f"no files for {key}")
    return total, edges


def slice_counts(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Counts in the 16 slices [5i, 5i + 5) of 0-80%; every slice edge must be a bin edge."""
    edges = np.asarray(edges, dtype=float)
    positions = np.searchsorted(edges, SLICE_EDGES)
    if np.any(positions >= len(edges)) or not np.array_equal(edges[positions], SLICE_EDGES):
        raise ValueError("the 5% slice edges of 0-80% must be bin edges of the centrality histogram")
    return np.add.reduceat(np.asarray(values, dtype=float)[positions[0]:positions[-1]], positions[:-1] - positions[0])


def trigger_histogram(paths, trigger: str) -> tuple[np.ndarray | None, np.ndarray | None, list]:
    """The centrality distribution of the data target from FillPhotonQA data outputs (one run each): the events of the
    named scaled trigger, on the fill's own centrality bins (None when no output has one). Returns the values, their
    edges and one record per output (its bit, whether the trigger was enabled, its counts per 5% slice)."""
    total, edges, records = None, None, []
    for path in paths:
        with uproot.open(path) as f:
            if 'meta/package' in f:
                metadata=json.loads(str(f['meta/package']))
                if metadata['system']!='AuAu' or metadata['family']!='data':
                    raise ValueError('centrality DATA target needs an AuAu DATA package')
                key='qa/events/centrality_scaled_bits'
                these=f[key].axes[0].edges()
                for source_index, raw in enumerate(f['records/sources']['json'].array(library='np')):
                    source=json.loads(raw); meta=source['metadata']
                    table=meta.get('trigger_selection') or meta.get('trigger_run_info')
                    if not table: raise ValueError('missing per-source trigger names in '+str(path))
                    bits=[b for b,name in enumerate(table['name']) if name==trigger]
                    if len(bits)>1: raise ValueError('ambiguous trigger name in package')
                    record=dict(path=str(path),source=source_index,input_sha256=source['sha256'],
                                bit=bits[0] if bits else None,enabled=bool(bits) and int(table['initial_prescale'][bits[0]])>=0)
                    if record['enabled']:
                        values=package.source_histograms(f,[key],source_index)[key][0][:,bits[0]]
                        if edges is None: total,edges=values.copy(),these
                        elif not np.array_equal(these,edges): raise ValueError('packaged centrality axes differ')
                        else: total+=values
                        record['counts']=slice_counts(values,these).tolist()
                    records.append(record)
                continue
            if "family" not in f or str(f["family"]) != "data" or "trigger_run_info" not in f:
                raise ValueError(f"{path} is not a FillPhotonQA data output")
            if str(f["trigger_run_info_consistent"]) != "1":
                raise ValueError(f"{path} mixes runs with different trigger tables; one run per output")
            table = json.loads(str(f["trigger_run_info"]))
            bits = [b for b, name in enumerate(table["name"]) if name == trigger]
            if len(bits) > 1:
                raise ValueError(f"{path}: the trigger {trigger!r} names more than one bit: {bits}")
            record = {"path": str(path), "sha256": sha256(path), "bit": bits[0] if bits else None,
                      "enabled": bool(bits) and int(table["initial_prescale"][bits[0]]) >= 0}
            if record["enabled"] and DATA_BITS in f:  # no histogram: no minimum-bias event in the vertex window
                h = f[DATA_BITS]
                values, these = np.asarray(h.values(), dtype=float)[:, bits[0]], h.axes[0].edges()
                if edges is None:
                    total, edges = values.copy(), these
                elif not np.array_equal(these, edges):
                    raise ValueError(f"{path}: {DATA_BITS} has other bin edges")
                else:
                    total += values
                record["counts"] = slice_counts(values, these).tolist()
            records.append(record)
    return total, edges, records


def trigger_counts(paths, trigger: str) -> tuple[np.ndarray, list]:
    """Data target counts per 5% slice (trigger_histogram) and the record of every output."""
    values, edges, records = trigger_histogram(paths, trigger)
    return (np.zeros(len(SLICE_EDGES) - 1) if values is None else slice_counts(values, edges)), records


def weights(data: np.ndarray, simulation: np.ndarray) -> np.ndarray:
    """p_DATA / p_SIM per slice, each normalised inside 0-80%; refuses an empty slice on either side (an empty data slice
    would give a zero weight, which FillPhotonJetResponse.C refuses: fail here, before any response job runs)."""
    data, simulation = np.asarray(data, dtype=float), np.asarray(simulation, dtype=float)
    if data.sum() <= 0 or simulation.sum() <= 0:
        raise ValueError("an empty centrality distribution")
    if np.any(simulation <= 0):
        raise ValueError(f"simulated slices without events: {np.flatnonzero(simulation <= 0).tolist()}")
    if np.any(data <= 0):
        raise ValueError(f"data slices without events: {np.flatnonzero(data <= 0).tolist()}")
    return (data / data.sum()) / (simulation / simulation.sum())


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], fromfile_prefix_chars="@")
    parser.add_argument("--config", required=True, type=Path, help="the configuration with the samples' normalisation")
    parser.add_argument("--data", required=True, nargs="+", type=Path, help="FillPhotonQA Au+Au data outputs, one run each")
    parser.add_argument("--data-trigger", default=DATA_TRIGGER, help="the scaled trigger of the data target, by its name in the trigger table")
    parser.add_argument("--simulation", required=True, nargs="+", action="append", metavar=("SAMPLE", "FILE"))
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.out.exists():
        raise FileExistsError(args.out)
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    data, runs = trigger_counts(args.data, args.data_trigger)
    simulation, per_sample = np.zeros(16), {}
    for group in args.simulation:
        name, paths = group[0], [Path(p) for p in group[1:]]
        weight = package.sample_weight(config["samples"][name], config.get("reference_cross_section_pb") or {})
        counts = slice_counts(*summed(paths, SIM_HISTOGRAM))
        per_sample[name] = {"weight_pb_per_event": weight, "counts": counts.tolist(), "files": [{"path": str(p), "sha256": sha256(p)} for p in paths]}
        simulation += weight * counts
    w = weights(data, simulation)
    result = {"schema": "AuAuCentralityWeightsV1", "formula": "w(c_i) = p_DATA(c_i) / p_SIM(c_i), 5% slices of 0-80%",
              "per_slice": w.tolist(), "argument": " ".join(repr(float(x)) for x in w),
              "data": {"target": "scaled trigger by name, minimum-bias classified, analysis vertex window", "histogram": DATA_BITS,
                       "trigger": args.data_trigger, "counts": data.tolist(), "runs_with_trigger": sum(r["enabled"] for r in runs),
                       "runs_without_trigger": sum(not r["enabled"] for r in runs), "files": runs},
              "simulation": {"histogram": SIM_HISTOGRAM, "counts": simulation.tolist(), "samples": per_sample}}
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    print(result["argument"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
