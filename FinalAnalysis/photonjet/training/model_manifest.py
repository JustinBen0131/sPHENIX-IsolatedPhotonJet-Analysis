"""Builder for a future complete model manifest; not wired into the V1 CLI."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

from photonjet.provenance import sha256_file


def build_model_manifest(
    *,
    name: str,
    system: str,
    model_path: Path,
    features: Iterable[str],
    score_direction: str,
    working_points: dict[str, Any],
    training_identity: str,
    split_identity: str,
) -> dict[str, Any]:
    feature_order = tuple(str(value) for value in features)
    if not name or system not in {"pp", "auau"}:
        raise ValueError("model name and supported system are required")
    if score_direction not in {"higher_is_signal", "lower_is_signal"}:
        raise ValueError("unsupported score direction")
    if not feature_order or len(set(feature_order)) != len(feature_order):
        raise ValueError("features must be nonempty and unique")
    model = Path(model_path).resolve()
    return {
        "schema": "PhotonJetModelManifestV1",
        "name": name,
        "system": system,
        "model": {
            "filename": model.name,
            "size_bytes": model.stat().st_size,
            "sha256": sha256_file(model),
        },
        "feature_order": list(feature_order),
        "score_direction": score_direction,
        "working_points": working_points,
        "training_identity": training_identity,
        "split_identity": split_identity,
    }
