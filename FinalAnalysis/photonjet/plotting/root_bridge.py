"""Generate a ROOT include from a verified plot contract."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from photonjet.plotting.contract import verify_plot_contract
from photonjet.provenance import sha256_file


def _cpp(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


def _header_text(contract: dict[str, Any]) -> str:
    contract = verify_plot_contract(contract)
    root = contract["annotations"]["root_tlatex"]
    lines = [root["dataset"], *root["cuts"]]
    if len(lines) > 7:
        raise ValueError("canonical ROOT header supports at most seven annotation lines")
    rendered = "\n".join(f'    "{_cpp(line)}",' for line in lines)
    return f'''// Generated from PhotonJetPlotContractV1. Do not edit physics labels.
#pragma once

#include <array>
#include <string_view>

namespace PhotonJetGeneratedAnnotation {{
struct Contract {{
inline static constexpr std::string_view kContractSha256 = "{contract["contract_sha256"]}";
inline static constexpr std::string_view kExperimentLabel = "{_cpp(root["experiment"])}";
inline static constexpr std::array<std::string_view, {len(lines)}> kLines = {{{{
{rendered}
}}}};
}};
inline constexpr Contract kAnnotation{{}};
}}  // namespace PhotonJetGeneratedAnnotation
'''


def emit_root_annotation_header(contract: dict[str, Any], output_path: Path) -> Path:
    text = _header_text(contract)
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")
    return output


def verify_root_annotation_header(
    contract: dict[str, Any], header_path: Path
) -> dict[str, Any]:
    """Regenerate in memory and verify that a generated include was not edited."""

    contract = verify_plot_contract(contract)
    header = Path(header_path).resolve()
    expected = _header_text(contract).encode("utf-8")
    if header.read_bytes() != expected:
        raise ValueError("generated ROOT annotation header was edited or is stale")
    return {
        "schema": "PhotonJetRootAnnotationReceiptV1",
        "contract_sha256": contract["contract_sha256"],
        "header_sha256": sha256_file(header),
        "verified": True,
    }
