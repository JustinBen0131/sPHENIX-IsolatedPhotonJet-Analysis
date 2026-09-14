"""Deterministic group-aware split assignment."""

from __future__ import annotations

import hashlib
from typing import Iterable


def _fraction(group_identity: str, seed: int) -> float:
    digest = hashlib.sha256(f"{seed}\x1f{group_identity}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") / float(1 << 64)


def assign_group_split(
    group_identities: Iterable[str],
    *,
    seed: int,
    train_fraction: float = 0.70,
    validation_fraction: float = 0.15,
) -> list[str]:
    if not 0 < train_fraction < 1:
        raise ValueError("train_fraction must be between zero and one")
    if not 0 <= validation_fraction < 1 - train_fraction:
        raise ValueError("validation_fraction leaves no held-out partition")
    train_edge = train_fraction
    validation_edge = train_fraction + validation_fraction
    result = []
    for identity in group_identities:
        value = _fraction(str(identity), int(seed))
        result.append(
            "train"
            if value < train_edge
            else "validation"
            if value < validation_edge
            else "holdout"
        )
    return result
