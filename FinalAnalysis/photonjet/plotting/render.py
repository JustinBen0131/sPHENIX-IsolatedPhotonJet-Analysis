"""Render provenance-compiled histogram contracts in canonical sPHENIX style."""

from __future__ import annotations

import json
from pathlib import Path
import platform
from typing import Any

import numpy as np

from photonjet.plotting.contract import verify_plot_contract
from photonjet.provenance import artifact, sha256_file, write_json


def _boxes_overlap(first: Any, second: Any, padding: float = 1.0) -> bool:
    return not (
        first.x1 + padding <= second.x0
        or second.x1 + padding <= first.x0
        or first.y1 + padding <= second.y0
        or second.y1 + padding <= first.y0
    )


def render_histogram(
    *,
    contract: dict[str, Any],
    histogram_path: Path,
    output_path: Path,
    receipt_path: Path,
) -> dict[str, Any]:
    """Render one marker-only xJgamma histogram and verify label geometry."""

    contract = verify_plot_contract(contract)
    histogram_file = Path(histogram_path).resolve()
    if sha256_file(histogram_file) != contract.get("histogram_sha256"):
        raise ValueError("histogram changed after plot-contract compilation")
    payload = json.loads(histogram_file.read_text(encoding="utf-8"))
    edges = np.asarray(payload["axis"]["edges"], dtype=float)
    sumw = np.asarray(payload["sumw"], dtype=float)
    sumw2 = np.asarray(payload["sumw2"], dtype=float)
    centers = 0.5 * (edges[1:] + edges[:-1])
    half_widths = 0.5 * np.diff(edges)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for plot rendering") from exc

    style = {
        "font.family": "DejaVu Sans",
        "font.style": "normal",
        "font.weight": "normal",
        "mathtext.fontset": "dejavusans",
        "mathtext.default": "it",
        "text.usetex": False,
        "font.size": 12,
        "axes.labelsize": 15,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.linewidth": 1.2,
        "xtick.top": True,
        "ytick.right": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
    with plt.rc_context(style):
        figure = plt.figure(figsize=(9.2, 6.6), dpi=160)
        axes = figure.add_axes([0.12, 0.12, 0.85, 0.85])
        artists: list[Any] = []
        annotations = contract["annotations"]["matplotlib"]
        artists.append(
            axes.text(
                0.025,
                0.975,
                annotations["experiment"],
                ha="left",
                va="top",
                fontsize=16,
                transform=axes.transAxes,
            )
        )
        artists.append(
            axes.text(
                0.025,
                0.905,
                annotations["dataset"],
                ha="left",
                va="top",
                fontsize=11.5,
                transform=axes.transAxes,
            )
        )
        for index, line in enumerate(annotations["cuts"]):
            artists.append(
                axes.text(
                    0.025,
                    0.835 - 0.065 * index,
                    line,
                    ha="left",
                    va="top",
                    fontsize=11.5,
                    transform=axes.transAxes,
                )
            )

        axes.errorbar(
            centers,
            sumw,
            xerr=half_widths,
            yerr=np.sqrt(sumw2),
            linestyle="none",
            marker="o",
            markersize=5.5,
            markerfacecolor="black",
            markeredgecolor="black",
            color="black",
            capsize=0,
        )
        axes.set_xlim(float(edges[0]), float(edges[-1]))
        upper_envelope = sumw + np.sqrt(sumw2)
        lower_envelope = sumw - np.sqrt(sumw2)
        data_top = float(np.max(upper_envelope)) if len(sumw) else 0.0
        data_bottom = min(0.0, float(np.min(lower_envelope))) if len(sumw) else 0.0
        data_span = data_top - data_bottom
        if data_span <= 0:
            data_span = 1.0
        # Start with an internal annotation band. The exact required headroom
        # is solved from rendered text boxes below; it is never a hard-coded
        # external canvas margin.
        axes.set_ylim(data_bottom, data_top + 0.67 * data_span)
        axes.set_xlabel(contract["observable"]["x_label_matplotlib"])
        axes.set_ylabel(contract["observable"]["y_label"])
        axes.grid(False)

        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        axes_box = axes.get_window_extent(renderer)
        artist_boxes = [artist.get_window_extent(renderer) for artist in artists]
        annotation_bottom_pixel = min(float(box.y0) for box in artist_boxes)
        desired_data_top_pixel = annotation_bottom_pixel - 8.0
        target_fraction = (
            (desired_data_top_pixel - float(axes_box.y0)) / float(axes_box.height)
        )
        if not 0.0 < target_fraction < 1.0:
            plt.close(figure)
            raise RuntimeError("annotations leave no readable in-frame data region")
        required_y_top = data_bottom + data_span / target_fraction
        if required_y_top > float(axes.get_ylim()[1]):
            axes.set_ylim(data_bottom, required_y_top)
            figure.canvas.draw()
            renderer = figure.canvas.get_renderer()
            axes_box = axes.get_window_extent(renderer)
            artist_boxes = [artist.get_window_extent(renderer) for artist in artists]
        inside = all(
            box.x0 >= axes_box.x0
            and box.y0 >= axes_box.y0
            and box.x1 <= axes_box.x1
            and box.y1 <= axes_box.y1
            for box in artist_boxes
        )
        data_top_pixel = float(axes.transData.transform((float(edges[0]), data_top))[1])
        annotation_bottom_pixel = min(float(box.y0) for box in artist_boxes)
        data_overlap = data_top_pixel + 6.0 > annotation_bottom_pixel
        annotation_overlap = any(
            _boxes_overlap(artist_boxes[left], artist_boxes[right])
            for left in range(len(artist_boxes))
            for right in range(left + 1, len(artist_boxes))
        )
        minimum_font = min(float(artist.get_fontsize()) for artist in artists)
        annotation_bottom_fraction = (
            (annotation_bottom_pixel - float(axes_box.y0)) / float(axes_box.height)
        )
        data_top_fraction = (
            (data_top_pixel - float(axes_box.y0)) / float(axes_box.height)
        )
        if not inside:
            plt.close(figure)
            raise RuntimeError("annotation geometry leaves the plotting canvas")
        if data_overlap:
            plt.close(figure)
            raise RuntimeError("annotation band overlaps the data/error envelope")
        if annotation_overlap:
            plt.close(figure)
            raise RuntimeError("annotation blocks overlap each other")
        if minimum_font < 11.0:
            plt.close(figure)
            raise RuntimeError("annotation font is below the audience-readable floor")

        output = Path(output_path).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=160, facecolor="white")
        plt.close(figure)

    receipt = {
        "schema": "PhotonJetPlotRenderReceiptV1",
        "artifact_state": "CANDIDATE_RAW",
        "contract_sha256": contract["contract_sha256"],
        "histogram_sha256": contract["histogram_sha256"],
        "image": artifact(output, "plot_image"),
        "rendered_annotations": contract["annotations"],
        "semantic_coverage": contract["semantic_coverage"],
        "semantic_coverage_sha256": contract["semantic_coverage_sha256"],
        "visible_semantic_keys": contract["visible_semantic_keys"],
        "renderer": {
            "implementation": "photonjet.plotting.render.render_histogram",
            "python": platform.python_version(),
            "matplotlib": matplotlib.__version__,
            "numpy": np.__version__,
            "backend": "Agg",
            "figure_size_inches": [9.2, 6.6],
            "dpi": 160,
        },
        "qa": {
            "all_annotations_inside_canvas": True,
            "all_annotations_inside_plotting_frame": True,
            "external_header_band": False,
            "annotation_data_overlap": False,
            "annotation_pair_overlap": False,
            "minimum_annotation_font_points": minimum_font,
            "marker_only_yield_rendering": True,
            "root_tlatex_experiment_label": contract["annotations"]["root_tlatex"]["experiment"],
            "annotation_band_bottom_axes_fraction": annotation_bottom_fraction,
            "data_envelope_top_axes_fraction": data_top_fraction,
            "annotation_data_padding_pixels": annotation_bottom_pixel - data_top_pixel,
        },
    }
    write_json(receipt_path, receipt)
    return receipt
