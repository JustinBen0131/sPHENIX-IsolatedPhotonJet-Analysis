"""Fail-closed promotion boundary for canonical plot artifacts.

Rendering produces a candidate.  Acceptance independently recompiles the
semantic contract from its histogram ancestors, replays the canonical
renderer, and requires exact image and receipt bytes before issuing the only
``ACCEPTED_CANONICAL`` token.  The public API has no free-form label or QA
override.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from photonjet.plotting.contract import compile_plot_contract
from photonjet.plotting.render import render_histogram
from photonjet.provenance import artifact, canonical_json, sha256_file, write_json


CANDIDATE_RAW = "CANDIDATE_RAW"
LEGACY_UNVERIFIED_ANNOTATION = "LEGACY_UNVERIFIED_ANNOTATION"
ACCEPTED_CANONICAL = "ACCEPTED_CANONICAL"
ACCEPTANCE_POLICY_ID = "canonical_plot_replay_and_on_canvas_semantics_v1"

_ACCEPTANCE_CHECKS = {
    "contract_recompiled_from_ancestors": True,
    "canonical_renderer_image_byte_exact": True,
    "canonical_renderer_receipt_byte_exact": True,
    "semantic_coverage_complete": True,
    "all_annotations_inside_plotting_frame": True,
    "external_header_band": False,
    "annotation_data_overlap": False,
    "annotation_pair_overlap": False,
    "minimum_annotation_font_points_at_least": 11.0,
}


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def _load_canonical_json(path: Path, description: str) -> dict[str, Any]:
    source = Path(path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    value = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{description} must be a JSON object")
    if source.read_bytes() != canonical_json(value).encode("utf-8"):
        raise ValueError(f"{description} is not canonical JSON")
    return value


def _audit_only_keys(contract: dict[str, Any]) -> list[str]:
    return sorted(
        str(entry["semantic_key"])
        for entry in contract["semantic_coverage"]["entries"]
        if entry["disposition"] == "audit_only"
    )


def _acceptance_body(
    *,
    contract: dict[str, Any],
    image_path: Path,
    render_receipt_path: Path,
) -> dict[str, Any]:
    return {
        "schema": "PhotonJetPlotAcceptanceReceiptV1",
        "state": ACCEPTED_CANONICAL,
        "source_state": CANDIDATE_RAW,
        "acceptance_policy_id": ACCEPTANCE_POLICY_ID,
        "contract_sha256": contract["contract_sha256"],
        "histogram_sha256": contract["histogram_sha256"],
        "histogram_receipt_sha256": contract["histogram_receipt_sha256"],
        "dataset_manifest_sha256": contract["dataset_manifest_sha256"],
        "style_id": contract["style_id"],
        "image": artifact(image_path, "accepted_plot_image"),
        "render_receipt": artifact(render_receipt_path, "plot_render_receipt"),
        "semantic_coverage_sha256": contract["semantic_coverage_sha256"],
        "visible_semantic_keys": contract["visible_semantic_keys"],
        "audit_only_semantic_keys": _audit_only_keys(contract),
        "checks": dict(_ACCEPTANCE_CHECKS),
    }


def _verify_canonical_replay(
    *,
    contract: dict[str, Any],
    histogram_path: Path,
    image_path: Path,
    render_receipt_path: Path,
) -> dict[str, Any]:
    supplied_render_receipt = _load_canonical_json(
        render_receipt_path, "plot render receipt"
    )
    if supplied_render_receipt.get("artifact_state") != CANDIDATE_RAW:
        raise ValueError("only a CANDIDATE_RAW canonical render can be accepted")

    # The replay uses the same basenames because artifact identities are
    # intentionally portable and path-free.  Exact equality then proves both
    # the pixels and every geometry/semantic receipt field came from the
    # current canonical renderer.
    with TemporaryDirectory(prefix="photonjet-plot-accept-") as directory:
        replay_root = Path(directory)
        replay_image = replay_root / image_path.name
        replay_receipt = replay_root / render_receipt_path.name
        render_histogram(
            contract=contract,
            histogram_path=histogram_path,
            output_path=replay_image,
            receipt_path=replay_receipt,
        )
        if image_path.read_bytes() != replay_image.read_bytes():
            raise ValueError("plot image differs from exact canonical renderer replay")
        if render_receipt_path.read_bytes() != replay_receipt.read_bytes():
            raise ValueError("plot render receipt differs from exact canonical renderer replay")

    qa = supplied_render_receipt.get("qa", {})
    if (
        qa.get("all_annotations_inside_canvas") is not True
        or qa.get("all_annotations_inside_plotting_frame") is not True
        or qa.get("external_header_band") is not False
        or qa.get("annotation_data_overlap") is not False
        or qa.get("annotation_pair_overlap") is not False
        or float(qa.get("minimum_annotation_font_points", 0.0)) < 11.0
        or float(qa.get("annotation_data_padding_pixels", 0.0)) < 6.0
    ):
        raise ValueError("plot render receipt does not satisfy canonical geometry policy")
    if supplied_render_receipt.get("semantic_coverage") != contract["semantic_coverage"]:
        raise ValueError("rendered semantic coverage differs from the active contract")
    return supplied_render_receipt


def _verify_acceptance_fields(
    *,
    acceptance: dict[str, Any],
    contract: dict[str, Any],
    image_path: Path,
    render_receipt_path: Path,
) -> dict[str, Any]:
    expected_keys = {
        "schema", "state", "source_state", "acceptance_policy_id",
        "contract_sha256", "histogram_sha256", "histogram_receipt_sha256",
        "dataset_manifest_sha256", "style_id", "image", "render_receipt",
        "semantic_coverage_sha256", "visible_semantic_keys",
        "audit_only_semantic_keys", "checks", "acceptance_sha256",
    }
    if set(acceptance) != expected_keys:
        missing = sorted(expected_keys - set(acceptance))
        extra = sorted(set(acceptance) - expected_keys)
        raise ValueError(f"plot acceptance fields differ; missing={missing}, extra={extra}")
    if acceptance["schema"] != "PhotonJetPlotAcceptanceReceiptV1":
        raise ValueError("unsupported plot acceptance schema")
    if acceptance["state"] != ACCEPTED_CANONICAL:
        raise ValueError("plot is not ACCEPTED_CANONICAL")
    if acceptance["source_state"] != CANDIDATE_RAW:
        raise ValueError("canonical acceptance has an invalid source state")
    if acceptance["acceptance_policy_id"] != ACCEPTANCE_POLICY_ID:
        raise ValueError("plot acceptance policy differs")
    body = {key: value for key, value in acceptance.items() if key != "acceptance_sha256"}
    if acceptance["acceptance_sha256"] != _sha256_json(body):
        raise ValueError("plot acceptance receipt changed after issuance")
    expected_body = _acceptance_body(
        contract=contract,
        image_path=image_path,
        render_receipt_path=render_receipt_path,
    )
    if body != expected_body:
        raise ValueError("plot acceptance differs from its recompiled semantic lineage")
    return acceptance


def accept_rendered_plot(
    *,
    histogram_path: Path,
    histogram_receipt_path: Path,
    dataset_manifest_path: Path,
    image_path: Path,
    render_receipt_path: Path,
    acceptance_receipt_path: Path,
) -> dict[str, Any]:
    """Replay and accept one exact canonical render.

    A bare PNG, a legacy annotation sidecar, a manually edited image, or a
    render made from stale ancestors fails before an acceptance token exists.
    """

    histogram = Path(histogram_path).resolve()
    histogram_receipt = Path(histogram_receipt_path).resolve()
    dataset_manifest = Path(dataset_manifest_path).resolve()
    image = Path(image_path).resolve()
    render_receipt = Path(render_receipt_path).resolve()
    output = Path(acceptance_receipt_path).resolve()
    inputs = {histogram, histogram_receipt, dataset_manifest, image, render_receipt}
    if output in inputs:
        raise ValueError("acceptance receipt must not overwrite an input artifact")
    if not image.is_file():
        raise FileNotFoundError(image)
    if not render_receipt.is_file():
        raise FileNotFoundError(render_receipt)

    contract = compile_plot_contract(
        histogram_path=histogram,
        histogram_receipt_path=histogram_receipt,
        dataset_manifest_path=dataset_manifest,
    )
    _verify_canonical_replay(
        contract=contract,
        histogram_path=histogram,
        image_path=image,
        render_receipt_path=render_receipt,
    )

    receipt = _acceptance_body(
        contract=contract,
        image_path=image,
        render_receipt_path=render_receipt,
    )
    receipt["acceptance_sha256"] = _sha256_json(receipt)
    write_json(output, receipt)
    acceptance = _load_canonical_json(output, "plot acceptance receipt")
    return _verify_acceptance_fields(
        acceptance=acceptance,
        contract=contract,
        image_path=image,
        render_receipt_path=render_receipt,
    )


def verify_plot_acceptance_receipt(
    *,
    acceptance_receipt_path: Path,
    histogram_path: Path,
    histogram_receipt_path: Path,
    dataset_manifest_path: Path,
    image_path: Path,
    render_receipt_path: Path,
) -> dict[str, Any]:
    """Rebuild and replay the complete portable plot-promotion proof."""

    acceptance = _load_canonical_json(
        Path(acceptance_receipt_path), "plot acceptance receipt"
    )
    histogram = Path(histogram_path).resolve()
    image = Path(image_path).resolve()
    render_receipt = Path(render_receipt_path).resolve()
    contract = compile_plot_contract(
        histogram_path=histogram,
        histogram_receipt_path=Path(histogram_receipt_path).resolve(),
        dataset_manifest_path=Path(dataset_manifest_path).resolve(),
    )
    _verify_canonical_replay(
        contract=contract,
        histogram_path=histogram,
        image_path=image,
        render_receipt_path=render_receipt,
    )
    return _verify_acceptance_fields(
        acceptance=acceptance,
        contract=contract,
        image_path=image,
        render_receipt_path=render_receipt,
    )
