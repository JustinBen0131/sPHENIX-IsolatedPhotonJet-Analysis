#!/usr/bin/env python3
"""Complete a bounded native SIM response grid; no event processing or dispatch.

Retain the caller's already admitted limits and existing recipe strings. The
native parser remains the authority for selection, model and resource checks.
"""
import argparse
import copy
import json
from pathlib import Path


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key: " + key)
        result[key] = value
    return result


def read(path):
    return json.loads(Path(path).read_text(), object_pairs_hook=unique_object)


def items(recipe):
    return unique_object(token.split("=", 1) for token in recipe.split())


def complete(base, contract, system, prescriptions=None):
    if system not in ("PP", "AuAu"):
        raise ValueError("system must be PP or AuAu")
    if set(base) != {"schema", "variants", "limits"} or base["schema"] != "PhotonJetNativeResponseGridSettingsV1":
        raise ValueError("base must be an unchanged ordinary native response grid")
    expected_limits = {"max_histogram_bytes", "max_resident_bytes", "max_output_bytes",
                       "max_sparse_bins_per_histogram", "max_seconds"}
    if set(base["limits"]) != expected_limits or any(type(v) is not int or v <= 0 for v in base["limits"].values()):
        raise ValueError("base requires all five explicit positive integer limits")
    if prescriptions is None:
        prescriptions = read(Path(__file__).parent / "config/response_grid_variations.json")
    expected = prescriptions[system]
    wanted = {v["name"]: v["variation"] for v in expected}
    found = set()
    for variant in base["variants"]:
        if set(variant) != {"name", "variation"}:
            raise ValueError("invalid base variant")
        name, recipe = variant["name"], variant["variation"]
        if name in found or name not in wanted or items(recipe) != items(wanted[name]):
            raise ValueError("duplicate, extra or changed base recipe: " + name)
        found.add(name)
    required = {"nominal", "smnom", "smjer_low", "smjer_high", "smemr_low", "smemr_high"}
    if not required <= found:
        raise ValueError("base must contain the complete existing six-component grid")
    if contract.get("schema") != "EmsNativeSelectionContractV1" or not {"source", "cone", "selection", "models", "pp_witness_line"} <= set(contract):
        raise ValueError("explicit source-bound EMS selection contract required")
    result = copy.deepcopy(base)
    result["variants"].extend(copy.deepcopy(v) for v in expected if v["name"] not in found)
    result["ems_selection_contract"] = copy.deepcopy(contract)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--ems-contract", type=Path, required=True)
    parser.add_argument("--system", choices=["PP", "AuAu"], required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = complete(read(args.base), read(args.ems_contract), args.system)
    with args.out.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"output": str(args.out), "components": [v["name"] for v in result["variants"]],
                      "execution_authorized": False}))


if __name__ == "__main__":
    main()
