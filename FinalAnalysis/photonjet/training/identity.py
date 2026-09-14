"""Stable, source-qualified identities for training rows and events."""

from __future__ import annotations

import hashlib
from typing import Any


def _digest(parts: tuple[Any, ...]) -> str:
    text = "\x1f".join(str(part) for part in parts)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def event_identity(source_sample: str, input_file_index: int, run: int, event: int) -> str:
    if not source_sample:
        raise ValueError("source_sample is required")
    return _digest(("event", source_sample, int(input_file_index), int(run), int(event)))


def candidate_identity(
    source_sample: str,
    input_file_index: int,
    run: int,
    event: int,
    input_tree_entry: int,
) -> str:
    return _digest(
        (
            "candidate",
            source_sample,
            int(input_file_index),
            int(run),
            int(event),
            int(input_tree_entry),
        )
    )
