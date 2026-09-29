#!/usr/bin/env python3
"""Compare the active environment with the pins in environment.yml.

    python TreeProduction/scripts/check_environment.py [environment.yml]

Prints one line per pinned package and exits nonzero when a package is missing or its version differs.
"""

from __future__ import annotations

import importlib.metadata
from pathlib import Path
import platform
import sys

import yaml


def pins(path: Path) -> dict[str, str]:
    """Package name -> pinned version, from the conda entries (name=version) and the pip entries (name==version)."""

    out: dict[str, str] = {}

    def add(spec: str, separator: str) -> None:
        name, _, version = spec.partition(separator)
        if version:
            out[name.strip().lower()] = version.strip()

    for item in yaml.safe_load(path.read_text(encoding="utf-8"))["dependencies"]:
        if isinstance(item, dict):
            for spec in item.get("pip", []):
                add(spec, "==")
        else:
            add(item, "=")
    return out


def installed(name: str) -> str | None:
    """The installed version of one pinned package, or None when it is missing."""

    if name == "python":
        return platform.python_version()
    if name == "root":
        try:
            import ROOT
        except ImportError:
            return None
        code = ROOT.gROOT.GetVersionInt()  # 63404 for 6.34.04
        return f"{code // 10000}.{code // 100 % 100}.{code % 100}"
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    path = Path(argv[0]) if argv else Path(__file__).resolve().parents[2] / 'environment.yml'
    failures = 0
    for name, pinned in pins(path).items():
        found = installed(name)
        ok = found == pinned
        failures += not ok
        print(f"{name:12s} {pinned:10s} {'ok' if ok else 'MISSING' if found is None else 'installed ' + found}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
