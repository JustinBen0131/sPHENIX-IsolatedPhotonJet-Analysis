"""Provenance-compiled plotting for PhotonJetTrees_v1 products."""

from photonjet.plotting.acceptance import (
    ACCEPTED_CANONICAL,
    CANDIDATE_RAW,
    LEGACY_UNVERIFIED_ANNOTATION,
    accept_rendered_plot,
    verify_plot_acceptance_receipt,
)
from photonjet.plotting.contract import compile_plot_contract, write_histogram_receipt
from photonjet.plotting.render import render_histogram

__all__ = [
    "ACCEPTED_CANONICAL",
    "CANDIDATE_RAW",
    "LEGACY_UNVERIFIED_ANNOTATION",
    "accept_rendered_plot",
    "compile_plot_contract",
    "render_histogram",
    "verify_plot_acceptance_receipt",
    "write_histogram_receipt",
]
