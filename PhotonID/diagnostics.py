"""Numerical payload first, rendering second. Rendering never trains or refits.

The families follow the recovered score/ROC, WP grid, centrality fit, and region
renderers. Numerical estimators live in numerics.py; layout is local and reusable.
Test/all panels evaluate frozen validation functions. ALL is descriptive only.
"""

from pathlib import Path
import numpy as np
from registry import write_json, read
from calibration import (
    PT_EDGES,
    CENT_EDGES,
    thresholds,
    efficiency_report,
    id_regions,
    line,
)


def save_diagnostics(directory, payload, wp, iso, population, model_name):
    from sklearn.metrics import roc_auc_score, roc_curve

    directory = Path(directory)
    directory.mkdir(parents=True)
    np.savez_compressed(directory / "payload.npz", **payload)
    score, y, w = payload["score"], payload["y"], payload["weight"]
    if set(y) != {0, 1}:
        raise ValueError(f"{population} diagnostics require both classes")
    if not np.all(np.isfinite(score)) or not np.all(np.isfinite(w) & (w > 0)):
        raise ValueError("invalid diagnostic scores/weights")
    fpr, tpr, _ = roc_curve(y, score, sample_weight=w)
    report = dict(
        population=population,
        label="full-sample / not held-out" if population == "all" else population,
        model_name=model_name,
        weighted_auc=float(roc_auc_score(y, score, sample_weight=w)),
        n=len(y),
        sumw=float(w.sum()),
        sumw2=float((w * w).sum()),
        roc=dict(fpr=fpr.tolist(), tpr=tpr.tolist()),
        achieved_efficiencies=efficiency_report(
            score, w, payload["et"], payload["centrality"], y == 1, wp
        ),
        wp=wp,
        isolation_package=iso,
    )
    if iso:
        report["isolation_achieved_efficiencies"] = efficiency_report(
            payload["isolation"],
            np.ones(len(payload["isolation"])),
            payload["isolation_et"],
            payload["isolation_centrality"],
            np.ones(len(payload["isolation"]), bool),
            iso,
        )
    if "reference_score" in payload:
        from scipy.stats import spearmanr

        ref = payload["reference_score"]
        if not np.all(np.isfinite(ref)):
            raise ValueError("nonfinite reference prediction")
        rf, rt, _ = roc_curve(y, ref, sample_weight=w)
        report["reference"] = dict(
            weighted_auc=float(roc_auc_score(y, ref, sample_weight=w)),
            pearson=float(np.corrcoef(ref, score)[0, 1]),
            spearman=float(spearmanr(ref, score).statistic),
            mean_absolute_delta=float(np.mean(np.abs(score - ref))),
            roc=dict(fpr=rf.tolist(), tpr=rt.tolist()),
        )
    write_json(directory / "numerical.json", report)
    render(directory)
    return {k: report[k] for k in ("weighted_auc", "n", "sumw", "sumw2", "label")}


def render(directory):
    """Regenerate every PNG using saved payload/tables, with zero model calls."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    directory = Path(directory)
    d = dict(np.load(directory / "payload.npz", allow_pickle=False))
    r = read(directory / "numerical.json")
    wp = r["wp"]
    auau = wp["system"] == "auau"
    score, y, w = d["score"], d["y"], d["weight"]
    title = f"{r['model_name']} — {r['label']} — candidate, not accepted"
    colors = ["#0072b2", "#d55e00"]

    def finish(fig, name):
        fig.suptitle(title, fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(directory / name, dpi=150)
        plt.close(fig)

    def hist(ax, values, weights, mask, label, bins, color):
        if np.any(mask) and weights[mask].sum() > 0:
            ax.hist(
                values[mask],
                bins=bins,
                weights=weights[mask] / weights[mask].sum(),
                histtype="step",
                label=label,
                color=color,
            )

    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for c, label in [(1, "Prompt photon signal"), (0, "Nonprompt background")]:
        hist(axs[0], score, w, y == c, label, np.linspace(0, 1, 51), colors[c])
    axs[0].set(xlabel="BDT score", ylabel="Training-recipe weighted fraction / bin")
    axs[0].legend()
    axs[1].plot(
        r["roc"]["fpr"], r["roc"]["tpr"], label=f"AUC = {r['weighted_auc']:.4f}"
    )
    axs[1].plot([0, 1], [0, 1], "k--", alpha=0.3)
    axs[1].set(xlabel="Background acceptance", ylabel="Signal efficiency")
    axs[1].legend()
    finish(fig, "01_score_distributions_and_roc.png")

    # Validation displays the actual quantile points and their uncertainties.
    # Other populations display achieved efficiencies at frozen curves.
    def grid(package, isolation, name):
        fig, axs = (
            plt.subplots(4, 4, figsize=(16, 12))
            if auau
            else plt.subplots(1, 1, figsize=(8, 5))
        )
        axes = np.atleast_1d(axs).ravel()
        val = d["isolation"] if isolation else score
        et = d["isolation_et"] if isolation else d["et"]
        cent = d["isolation_centrality"] if isolation else d["centrality"]
        weights = np.ones(len(val)) if isolation else w
        signal = np.ones(len(val), bool) if isolation else y == 1
        for i, ax in enumerate(axes):
            lo, hi = (5 * i, 5 * i + 5) if auau else (0, 1)
            for target in (70, 80, 90):
                key = ("ISO" if isolation else "WP") + str(target)
                if r["population"] == "validation":
                    cells = [
                        v
                        for v in package["tables"]["cells"]
                        if round(100 * v["efficiency"]) == target
                        and (not auau or v["cent_min"] == lo)
                    ]
                    ax.errorbar(
                        [(v["pt_min"] + v["pt_max"]) / 2 for v in cells],
                        [v["threshold"] for v in cells],
                        yerr=[v["error"] for v in cells],
                        marker="o",
                        label=key,
                    )
                    xx = np.linspace(15, 35, 100)
                    ax.plot(
                        xx,
                        line(
                            package["curves"][key],
                            np.full(100, (lo + hi) / 2) if auau else xx,
                        ),
                        "--",
                        alpha=0.5,
                    )
                else:
                    cut = thresholds(package, et, cent)[key]
                    eff = []
                    for a, b in zip(PT_EDGES[:-1], PT_EDGES[1:]):
                        m = (
                            signal
                            & np.isfinite(val)
                            & np.isfinite(cut)
                            & (et >= a)
                            & (et < b)
                        )
                        if auau:
                            m &= (cent >= lo) & (cent < hi)
                        passed = val < cut if isolation else val > cut
                        eff.append(
                            float(weights[m & passed].sum() / weights[m].sum())
                            if weights[m].sum() > 0
                            else np.nan
                        )
                    ax.plot(
                        [(a + b) / 2 for a, b in zip(PT_EDGES[:-1], PT_EDGES[1:])],
                        eff,
                        "o-",
                        label=key,
                    )
                    ax.axhline(target / 100, ls=":", alpha=0.3)
                ax.set(
                    xlabel="Photon ET [GeV]",
                    ylabel=(
                        ("Isolation threshold [GeV]" if isolation else "BDT threshold")
                        if r["population"] == "validation"
                        else "Efficiency at frozen cut"
                    ),
                )
            if auau:
                ax.set_title(f"{lo}–{hi}%")
            if i == 0:
                ax.legend(fontsize=8)
        finish(fig, name)

    def fits(package, isolation, name):
        fig, ax = plt.subplots(figsize=(9, 5))
        x = np.linspace(*package["domain"], 200)
        for f in package["tables"]["linear_fits"]:
            ax.plot(
                x,
                line(package["curves"][f["name"]], x),
                label=f["name"] + " frozen validation fit",
            )
            rows = f["inputs"]
            ax.errorbar(
                [v["axis"] for v in rows],
                [v["threshold"] for v in rows],
                yerr=[v["error"] for v in rows],
                fmt="o",
                alpha=0.5,
            )
        ax.set(
            xlabel="Centrality [%]" if auau else "Photon ET [GeV]",
            ylabel="Isolation threshold [GeV]" if isolation else "BDT threshold",
        )
        ax.legend()
        if r["population"] != "validation":
            ax.set_title(
                "Validation calibration reference; no refit on this population"
            )
        finish(fig, name)

    grid(
        wp,
        False,
        "02_bdt_wp_pt_centrality_grid.png" if auau else "02_bdt_wp_pt_grid.png",
    )
    fits(
        wp, False, "03_bdt_wp_centrality_fits.png" if auau else "03_bdt_wp_pt_fits.png"
    )
    fig, axs = plt.subplots(1, 3 if auau else 1, figsize=(14 if auau else 8, 4))
    cuts = thresholds(wp, d["et"], d["centrality"])
    regions = id_regions(score, cuts["WP70"], cuts["WP80"])
    for ax, (lo, hi) in zip(
        np.atleast_1d(axs), [(0, 20), (20, 50), (50, 80)] if auau else [(15, 35)]
    ):
        axis = d["centrality"] if auau else d["et"]
        m = (axis >= lo) & (axis < hi)
        for c, label in [(1, "Prompt"), (0, "Nonprompt")]:
            hist(ax, score, w, m & (y == c), label, np.linspace(0, 1, 51), colors[c])
        x = (lo + hi) / 2
        upper = line(wp["curves"]["WP80"], x)
        tight = line(wp["curves"]["WP70"], x)
        ax.axvspan(0.1, upper, color="#e69f00", alpha=0.12)
        ax.axvspan(upper, tight, color="grey", alpha=0.15)
        ax.axvspan(tight, 1, color="#0072b2", alpha=0.1)
        counts = [int(np.sum(m & (regions == v))) for v in (2, 3, 4)]
        ax.set(
            title=f"{lo}–{hi}"
            + ("%" if auau else " GeV")
            + f"\nNT / gap / tight: {counts}",
            xlabel="BDT score",
            ylabel="Weighted class fraction / bin",
        )
        ax.legend(fontsize=8)
    fig.text(
        0.5,
        0.01,
        "Shading uses bin midpoint; counts use each candidate’s frozen threshold. No reconstructed-isolation cut.",
        ha="center",
        fontsize=8,
    )
    finish(fig, "04_bdt_score_regions.png")
    if r["isolation_package"]:
        grid(r["isolation_package"], True, "05_isolation_wp_pt_centrality_grid.png")
        fits(r["isolation_package"], True, "06_isolation_wp_centrality_fits.png")
    if "reference" in r:
        ref = d["reference_score"]
        rr = r["reference"]
        fig, axs = plt.subplots(1, 3, figsize=(15, 4))
        axs[0].plot(
            rr["roc"]["fpr"],
            rr["roc"]["tpr"],
            label=f"Original AUC {rr['weighted_auc']:.4f}",
        )
        axs[0].plot(
            r["roc"]["fpr"], r["roc"]["tpr"], label=f"Retrained {r['weighted_auc']:.4f}"
        )
        axs[0].legend()
        axs[0].set(xlabel="Background acceptance", ylabel="Signal efficiency")
        axs[1].hist2d(
            ref,
            score,
            bins=80,
            range=[[0, 1], [0, 1]],
            norm=matplotlib.colors.LogNorm(),
        )
        axs[1].set(
            xlabel="Original PPG12 score",
            ylabel="Retrained score",
            title=f"Pearson {rr['pearson']:.4f}; Spearman {rr['spearman']:.4f}",
        )
        for c, label in [(1, "Prompt"), (0, "Nonprompt")]:
            hist(
                axs[2],
                score - ref,
                w,
                y == c,
                label,
                np.linspace(-1, 1, 101),
                colors[c],
            )
        axs[2].set(
            xlabel="Retrained − original score", ylabel="Weighted class fraction / bin"
        )
        axs[2].legend()
        finish(fig, "07_ppg12_reference_comparison.png")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Render one saved PhotonID diagnostics population; no training or fitting."
    )
    p.add_argument("directory", type=Path)
    render(p.parse_args().directory)
