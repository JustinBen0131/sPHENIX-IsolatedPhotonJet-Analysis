"""Small, deterministic provenance records shared by every public stage."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any, Iterable


def canonical_json(value: Any) -> str:
    """Serialize a numeric payload deterministically.

    JSON payloads never contain NaN or infinity. Those values are analysis
    failures, not portable data.
    """

    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ) + "\n"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit(repo: Path) -> str | None:
    completed = subprocess.run(
        ["git", "-C", str(Path(repo).resolve()), "rev-parse", "HEAD"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def artifact(path: Path, role: str) -> dict[str, Any]:
    """Bind portable artifact identity without serializing a machine location."""

    supplied = Path(path)
    resolved = supplied.resolve()
    return {
        "role": role,
        "name": supplied.name,
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def run_manifest(
    *,
    stage: str,
    inputs: Iterable[Path],
    outputs: Iterable[Path],
    parameters: dict[str, Any],
    validators: Iterable[str] = (),
    repo: Path | None = None,
) -> dict[str, Any]:
    return {
        "schema": "PhotonJetRunManifestV1",
        "stage": stage,
        "software_commit": git_commit(repo) if repo else None,
        "parameters": parameters,
        "inputs": [artifact(path, "input") for path in inputs],
        "outputs": [artifact(path, "output") for path in outputs],
        "validators": sorted(set(validators)),
    }


def write_json(path: Path, value: Any) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(canonical_json(value), encoding="utf-8")
