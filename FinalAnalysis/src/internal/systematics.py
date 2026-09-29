#!/usr/bin/env python3
"""Systematic uncertainties of the per-photon result, from measurement.py results of one change each.

    ./FinalAnalysis/final_analysis.sh runSystematics --config systematics.yaml --system auau --out results/auau_systematics.json
    ./FinalAnalysis/final_analysis.sh runSystematics --config systematics.yaml --ratio auau pp --out results/ratio_systematics.json

systematics.yaml names the nominal result of each system and its sources:

    nominal: {auau: results/auau.json, pp: results/pp.json}
    sources:
      jet_energy_scale:          # a two-sided nuisance: up and down templates
        family: jet_response
        correlation: jes_shared  # the same physical input in both systems: varied together in the ratio
        auau: {up: results/auau_jes_up.json, down: results/auau_jes_down.json}
        pp: {up: results/pp_jes_up.json, down: results/pp_jes_down.json}
      prior_form:                # a one-sided alternative, mirrored
        family: unfolding
        auau: {alternative: results/auau_prior_exp2.json}
      iterations:                # alternatives that probe one effect: the envelope is the display band
        family: unfolding
        auau: {alternatives: {k3: results/auau_k3.json, k4: results/auau_k4.json}, representative: k4}

Every source keeps its signed shift templates (varied minus nominal, per reported xJ bin). A two-sided source enters the
covariance as d d^T with d = (up - down) / 2, and the asymmetric bands take, per bin, the larger positive and the larger
negative of its two shifts. A one-sided source is mirrored: d = alternative - nominal, both signs. An envelope source shows the
largest |shift| over its alternatives as its band; its covariance is d d^T of one named representative (an envelope does not
define bin-to-bin correlations). Sources are independent of each other and add in covariance and in quadrature per band. The
statistical covariance is the nominal result's own (data toys plus response toys), counted once.

Ratio (system A / system B, per reported bin): a source present in both systems with the same correlation tag in both is
varied coherently, R(A + d_A, B + d_B) - R(A, B); every other source of either system varies its own system alone. A tag is
one string for both systems or one per system ({pp: tag, auau: tag}). Statistics of the two systems are independent.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import yaml

KINDS = ("updown", "one_sided", "envelope")


@dataclass(frozen=True)
class Result:
    """The reported part of one measurement.py result."""

    central: np.ndarray            # per-photon value per reported xJ bin
    statistical: np.ndarray        # its statistical covariance (data + response toys)
    xj_edges: tuple[float, ...]    # of the reported bins
    source: str

    @classmethod
    def from_json(cls, path: Path) -> "Result":
        doc = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_document(doc, str(path))

    @classmethod
    def from_document(cls, doc: Mapping[str, Any], source: str = "") -> "Result":
        shown = np.asarray(doc["reported_xj_mask"], dtype=bool)
        entry = doc["entry"]
        covariance = (np.asarray(entry["covariance_data"], dtype=float) + np.asarray(entry["covariance_response"], dtype=float))
        if covariance.shape[0] != shown.sum():  # stored on every bin
            covariance = covariance[np.ix_(shown, shown)]
        edges = np.asarray(doc["settings"]["xj_edges"], dtype=float)
        reported_edges = tuple(float(x) for x in np.r_[edges[:-1][shown], edges[1:][shown][-1:]])
        return cls(np.asarray(doc["per_photon"], dtype=float)[shown], covariance, reported_edges, source)


def _same_binning(nominal: Result, other: Result) -> None:
    if nominal.central.shape != other.central.shape or not np.allclose(nominal.xj_edges, other.xj_edges, rtol=0, atol=1e-12):
        raise ValueError(f"{other.source}: reported bins differ from the nominal's")


def source_shifts(nominal: Result, spec: Mapping[str, Any], load) -> dict[str, Any]:
    """Signed templates, covariance template and bands of one source of one system. ``load(path)`` returns a Result."""

    if "up" in spec or "down" in spec:
        up, down = load(spec["up"]), load(spec["down"])
        for r in (up, down):
            _same_binning(nominal, r)
        d_up, d_down = up.central - nominal.central, down.central - nominal.central
        return {"kind": "updown", "templates": {"up": d_up, "down": d_down}, "covariance_template": 0.5 * (d_up - d_down),
                "band_up": np.maximum(np.maximum(d_up, d_down), 0.0), "band_down": np.minimum(np.minimum(d_up, d_down), 0.0)}
    if "alternative" in spec:
        alternative = load(spec["alternative"])
        _same_binning(nominal, alternative)
        d = alternative.central - nominal.central
        return {"kind": "one_sided", "templates": {"alternative": d}, "covariance_template": d,
                "band_up": np.abs(d), "band_down": -np.abs(d)}
    if "alternatives" in spec:
        templates = {}
        for name, path in spec["alternatives"].items():
            r = load(path)
            _same_binning(nominal, r)
            templates[name] = r.central - nominal.central
        representative = spec.get("representative")
        if representative is None or representative not in templates:
            raise ValueError("an envelope source names one of its alternatives as the representative of its covariance")
        width = np.max(np.abs(np.vstack(list(templates.values()))), axis=0)
        return {"kind": "envelope", "templates": templates, "covariance_template": templates[representative],
                "representative": representative, "band_up": width, "band_down": -width}
    raise ValueError("a source gives up and down, an alternative, or alternatives")


def combine(nominal: Result, sources: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """Total systematic covariance, asymmetric bands and the statistical covariance of one system."""

    n = nominal.central.size
    covariance = np.zeros((n, n))
    up, down = np.zeros(n), np.zeros(n)
    for s in sources.values():
        d = np.asarray(s["covariance_template"], dtype=float)
        covariance += np.outer(d, d)
        up += np.square(s["band_up"])
        down += np.square(s["band_down"])
    statistical = nominal.statistical
    total = covariance + statistical
    sigma = np.sqrt(np.diag(total))
    return {"central": nominal.central, "xj_edges": nominal.xj_edges,
            "statistical_sigma": np.sqrt(np.diag(statistical)), "systematic_sigma": np.sqrt(np.diag(covariance)),
            "systematic_band_up": np.sqrt(up), "systematic_band_down": -np.sqrt(down),
            "systematic_covariance": covariance, "statistical_covariance": statistical,
            "total_correlation": total / np.outer(sigma, sigma) if np.all(sigma > 0) else None}


def correlation_tag(source: Mapping[str, Any], system: str):
    """A source's correlation tag in one system: one string for both systems, or one per system; None: system-specific."""

    tag = source.get("correlation")
    return tag.get(system) if isinstance(tag, Mapping) else tag


def ratio(a: Result, b: Result, shifts_a: Mapping[str, Mapping[str, Any]], shifts_b: Mapping[str, Mapping[str, Any]],
          tags: Mapping[str, Any]) -> dict[str, Any]:
    """Systematics of the per-bin ratio a / b. ``tags`` maps a source name to its correlation tag, shared by both systems,
    or to a (tag in a, tag in b) pair; the source varies coherently only when both tags are the same and not None."""

    _same_binning(a, b)
    central = a.central / b.central
    grouped: dict[str, dict[str, Any]] = {}
    for name in sorted(set(shifts_a) | set(shifts_b)):
        tag = tags.get(name)
        tag_a, tag_b = tag if isinstance(tag, tuple) else (tag, tag)
        coherent = tag_a is not None and tag_a == tag_b and name in shifts_a and name in shifts_b
        for side, shifts in (("a", shifts_a), ("b", shifts_b)):
            if name not in shifts:
                continue
            key = f"{name}" if coherent else f"{name}[{side}]"
            grouped.setdefault(key, {"a": np.zeros_like(central), "b": np.zeros_like(central), "coherent": coherent})
            grouped[key][side] = np.asarray(shifts[name]["covariance_template"], dtype=float)
    templates = {key: (a.central + g["a"]) / (b.central + g["b"]) - central for key, g in grouped.items()}
    covariance = sum((np.outer(d, d) for d in templates.values()), np.zeros((central.size, central.size)))
    # independent statistics, linearized: J_a V_a J_a^T + J_b V_b J_b^T
    ja, jb = np.diag(1.0 / b.central), np.diag(-a.central / b.central ** 2)
    statistical = ja @ a.statistical @ ja.T + jb @ b.statistical @ jb.T
    return {"central": central, "xj_edges": a.xj_edges, "templates": templates, "coherent": sorted(k for k, g in grouped.items() if g["coherent"]),
            "systematic_covariance": covariance, "systematic_sigma": np.sqrt(np.diag(covariance)),
            "statistical_covariance": statistical, "statistical_sigma": np.sqrt(np.diag(statistical))}


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _jsonable(value):
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def system_sources(config: Mapping[str, Any], system: str, nominal: Result, load) -> dict[str, dict[str, Any]]:
    out = {}
    for name, source in (config.get("sources") or {}).items():
        if system in source:
            out[name] = {**source_shifts(nominal, source[system], load), "family": source.get("family"),
                         "correlation": correlation_tag(source, system)}
    return out


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, type=Path)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--system")
    group.add_argument("--ratio", nargs=2, metavar=("NUMERATOR", "DENOMINATOR"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise FileExistsError(args.output)
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    base = args.config.parent
    used: dict[str, str] = {}

    def load(path) -> Result:
        full = (base / path) if not Path(path).is_absolute() else Path(path)
        used[str(full)] = _sha256(full)
        return Result.from_json(full)

    if args.system:
        nominal = load(config["nominal"][args.system])
        sources = system_sources(config, args.system, nominal, load)
        output = {"schema": "PhotonJetSystematicsV1", "system": args.system, **combine(nominal, sources), "sources": sources}
    else:
        a_name, b_name = args.ratio
        a, b = load(config["nominal"][a_name]), load(config["nominal"][b_name])
        shifts_a, shifts_b = system_sources(config, a_name, a, load), system_sources(config, b_name, b, load)
        tags = {name: (correlation_tag(source, a_name), correlation_tag(source, b_name))
                for name, source in (config.get("sources") or {}).items()}
        output = {"schema": "PhotonJetRatioSystematicsV1", "ratio": [a_name, b_name], **ratio(a, b, shifts_a, shifts_b, tags)}
    output["config"] = {"path": str(args.config), "sha256": _sha256(args.config)}
    output["inputs"] = used
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(_jsonable(output), indent=1, allow_nan=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
