"""
Build PSB 2027 figures F1-F4 from saved analysis outputs.

Design constraints (12-page limit, four figures):
  * each figure must fit ~0.45-0.55 page => single-column width 4.6in,
    height 2.0-2.6in
  * greyscale-safe: distinguish by shape/hatch/position, not colour alone
  * no chartjunk; the null band in F4 must be unambiguous at a glance

Outputs: figures/{F1_designs,F2_stability,F3_esm2,F4_nullband_*}.pdf
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Repository-relative paths. Scripts live in analysis/; everything they read
# and write is inside this repository.
REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUT = REPO / "outputs"
CACHE = REPO / "cache"
for _d in (OUT, CACHE):
    _d.mkdir(parents=True, exist_ok=True)
BASE = REPO  # legacy alias



# Deterministic PDF output: matplotlib stamps a CreationDate into every PDF,
# so otherwise-identical figures differ byte-for-byte on each run and any
# reproducibility check comparing committed artefacts fails spuriously.
DETERMINISTIC_PDF = {"CreationDate": None}

plt.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 8.5,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

W = 4.6  # single-column width, inches
GREY, DARK, ACC = "0.72", "0.25", "0.0"


# ---------------------------------------------------------------- F2
def fig2_stability():
    """Caller robustness: containment / rho / top-50, viable vs degenerate."""
    df = pd.read_csv(OUT / "E3_calibrated_summary.csv")
    models = ["Geneformer", "scGPT", "scFoundation"]
    viable = ["MAD z>3", "MAD z>3.5", "Top-n by MAD score"]
    degen = ["IQR k=3 (Tukey extreme)"]

    fig, axes = plt.subplots(1, 3, figsize=(W, 1.95), sharey=True)
    for ax, m in zip(axes, models):
        sub = df[df.model == m].set_index("method")
        methods = viable + degen
        n_orig = sub.loc["|z|>3 (original)", "n_outliers"]
        vals = [sub.loc[k, "containment"] if k in sub.index else np.nan
                for k in methods]
        # Containment is capped at min(n_new, n_orig)/n_orig: a stricter caller
        # returning fewer genes CANNOT contain them all. Plot the cap so
        # saturation is not misread as instability.
        caps = [min(sub.loc[k, "n_outliers"], n_orig) / n_orig
                if k in sub.index else np.nan for k in methods]
        xs = np.arange(len(methods))
        colors = [DARK] * len(viable) + [GREY] * len(degen)
        hatches = [""] * len(viable) + ["///"] * len(degen)
        for x, v, c, h in zip(xs, vals, colors, hatches):
            ax.bar(x, v, color=c, hatch=h, edgecolor="black", linewidth=0.5,
                   width=0.72)
        ax.scatter(xs, caps, marker="_", s=90, color="black", linewidth=1.1,
                   zorder=5)
        ax.set_xticks(xs)
        ax.set_xticklabels(["MAD\n>3", "MAD\n>3.5", "rank\nMAD", "IQR\n(deg.)"],
                           fontsize=6.5)
        ax.set_ylim(0, 1.08)
        ax.axhline(1.0, color="black", lw=0.5, ls=":")
        rho = sub["spearman_rho"].iloc[0]
        ax.set_title(f"{m}\n$\\rho$ = {rho:.3f}", fontsize=7.5)
    axes[0].set_ylabel("containment of\noriginal outliers")
    fig.tight_layout(pad=0.3)
    fig.savefig(REPO / "figures" / "F2_stability.pdf", metadata=DETERMINISTIC_PDF)
    plt.close(fig)
    print("  F2_stability.pdf")


# ---------------------------------------------------------------- F3
def fig3_esm2():
    """scFM-only fraction per class; mitochondrial is the informative exception."""
    d = pd.read_csv(OUT / "E6_scfm_only_by_class.csv")
    order = ["constrained", "disease", "ribosomal", "mitochondrial"]
    d = d.set_index("cls").loc[order].reset_index()
    frac = d.n_scfm_only / d.n_scfm_outliers

    fig, ax = plt.subplots(figsize=(W, 1.55))
    ys = np.arange(len(d))[::-1]
    for y, f, row in zip(ys, frac, d.itertuples()):
        is_mito = row.cls == "mitochondrial"
        ax.barh(y, f, color=GREY if is_mito else DARK,
                hatch="///" if is_mito else "", edgecolor="black",
                linewidth=0.5, height=0.62)
        ax.text(f + 0.015, y, f"{row.n_scfm_only}/{row.n_scfm_outliers}",
                va="center", fontsize=7)
    ax.set_yticks(ys)
    ax.set_yticklabels(["constrained", "disease\n(ClinVar)", "ribosomal",
                        "mitochondrial"])
    ax.set_xlim(0, 1.16)
    ax.set_xlabel("fraction of Geneformer outliers that are NOT ESM-2 outliers")
    ax.axvline(1.0, color="black", lw=0.5, ls=":")
    fig.tight_layout(pad=0.3)
    fig.savefig(REPO / "figures" / "F3_esm2.pdf", metadata=DETERMINISTIC_PDF)
    plt.close(fig)
    print("  F3_esm2.pdf")


# ---------------------------------------------------------------- F4
def fig4_nullband():
    """THE figure. Treatment delta vs matched-control null. macro-F1 only.

    One float, three panels: PBMC3k primary, PBMC3k sensitivity, and the
    Tabula Sapiens primary arm. The axes are deliberately NOT shared. The two
    datasets have different baselines (0.924 against 0.635) and different
    control counts (200 against 100), so a common x-axis would invite reading
    the three histograms as one null distribution. Each panel therefore states
    its own baseline and n.

    Values are plotted in units of 1e-3 so the tick labels fit a narrow panel.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    # dataset, control key, treatment key, title, n genes, note
    PANELS = [
        ("pbmc3k", "control_results_full", "treatment",
         "PBMC3k\nprimary", 50, None),
        ("pbmc3k", "control_results_no_ribo_mito", "sensitivity",
         "PBMC3k\nsensitivity", 36, "no ribo/mito"),
        ("tabula_sapiens", "control_results_full", "treatment",
         "Tabula Sapiens\nprimary only", 50, None),
    ]

    cache = {}

    def load(ds):
        if ds not in cache:
            path = OUT / f"E2_ablation_{ds}.json"
            if not path.exists():
                return None
            res = json.load(open(path))
            bp = OUT / f"E2_baseline_{ds}.json"
            base = (json.load(open(bp))["baseline_retrained_f1"]
                    if bp.exists() else res["baseline"]["baseline_retrained_f1"])
            cache[ds] = (res, base)
        return cache[ds]

    for ds, *_ in PANELS:
        if load(ds) is None:
            print(f"  (skipping F4: no E2_ablation_{ds}.json)")
            return

    fig, axes = plt.subplots(1, 3, figsize=(W, 2.05))
    for ax, (ds, ck, tk, title, ngenes, note) in zip(axes, PANELS):
        res, base = load(ds)
        ctrls = res[ck]
        cd = (np.array([c["retrained_f1"] for c in ctrls]) - base) * 1e3
        td = (res[tk]["retrained_f1"] - base) * 1e3
        lo, hi = np.percentile(cd, [2.5, 97.5])
        z = (td - cd.mean()) / cd.std(ddof=1)

        ax.axvspan(lo, hi, color="0.88", zorder=0)
        ax.hist(cd, bins=22 if len(ctrls) >= 200 else 16, color=GREY,
                edgecolor="black", linewidth=0.3, zorder=2)
        ax.axvline(td, color=ACC, lw=1.5, zorder=4)
        for b in (lo, hi):
            ax.axvline(b, color=DARK, lw=0.7, ls="--", zorder=3)

        ax.set_ylim(0, ax.get_ylim()[1] * 1.72)   # headroom for the annotation

        lines = [f"$z$ = {z:+.2f}  INSIDE"]
        if note:
            lines.append(note)
        lines += [f"{ngenes} genes, $n$ = {len(ctrls)}",
                  f"baseline {base:.4f}"]
        ax.text(0.05, 0.98, "\n".join(lines), transform=ax.transAxes,
                fontsize=5.5, va="top", linespacing=1.4, zorder=6,
                bbox=dict(boxstyle="round,pad=0.24", fc="white", ec="0.6",
                          lw=0.45))

        ax.set_title(title, fontsize=7.2)
        ax.tick_params(axis="x", pad=1.5)
        ax.locator_params(axis="x", nbins=4)

    axes[0].set_ylabel("matched control sets", fontsize=7.2)
    axes[1].set_xlabel("$\\Delta$ macro-$F_1$ vs baseline ($\\times 10^{-3}$)")

    fig.legend(handles=[Line2D([], [], color=ACC, lw=1.5, label="treatment"),
                        Patch(facecolor="0.88", label="central 95% of null"),
                        Line2D([], [], color=DARK, lw=0.7, ls="--",
                               label="2.5th / 97.5th percentile")],
               loc="lower center", ncol=3, frameon=False, fontsize=6,
               bbox_to_anchor=(0.5, -0.05), handlelength=1.5,
               columnspacing=1.3, handletextpad=0.5)
    fig.tight_layout(pad=0.3, w_pad=0.75, rect=(0, 0.05, 1, 1))
    fig.savefig(REPO / "figures" / "F4_nullband.pdf",
                metadata=DETERMINISTIC_PDF)
    plt.close(fig)
    print("  F4_nullband.pdf")


# ---------------------------------------------------------------- F1
def fig1_designs():
    """Three model designs, the shared geometry screen, and follow-up tests."""
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    width, height = 5.6, 2.85
    grey, edge = "0.35", "0.25"
    fills = ["#EAEFF5", "#EAF2EA", "#F5EEE6"]
    models = [
        {
            "name": "Geneformer",
            "design": "BERT-style encoder",
            "detail": [
                "rank-ordered gene tokens",
                "no expression values",
                "20,271 genes",
            ],
            "out": "410 outliers",
        },
        {
            "name": "scGPT",
            "design": "generative transformer",
            "detail": [
                "gene tokens $+$ expression",
                "value embeddings",
                "60,694 genes",
            ],
            "out": "188 outliers",
        },
        {
            "name": "scFoundation",
            "design": "asymmetric encoder-decoder",
            "detail": [
                "continuous expression via",
                "learned auto-discretisation",
                "19,264 genes",
            ],
            "out": "164 outliers",
        },
    ]

    def box(ax, x, y, w, h, fc, lw=0.7, rounding=0.015):
        ax.add_patch(FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle=f"round,pad=0,rounding_size={rounding}",
            linewidth=lw,
            edgecolor=edge,
            facecolor=fc,
            zorder=2,
        ))

    def arrow(ax, x0, y0, x1, y1, lw=0.7):
        ax.add_patch(FancyArrowPatch(
            (x0, y0),
            (x1, y1),
            arrowstyle="-|>",
            mutation_scale=6,
            linewidth=lw,
            color=edge,
            zorder=3,
            shrinkA=0,
            shrinkB=0,
        ))

    with plt.rc_context({
        "font.size": 7,
        "font.family": "serif",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.6,
        "savefig.bbox": None,
    }):
        fig, ax = plt.subplots(figsize=(width, height))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

        # Three model designs.
        bw, bh, gap = 0.30, 0.34, 0.05
        y0 = 0.60
        for i, model in enumerate(models):
            x = i * (bw + gap)
            box(ax, x, y0, bw, bh, fills[i])
            ax.text(
                x + bw / 2,
                y0 + bh - 0.055,
                model["name"],
                ha="center",
                va="top",
                fontsize=8,
                fontweight="bold",
                zorder=4,
            )
            ax.text(
                x + bw / 2,
                y0 + bh - 0.125,
                model["design"],
                ha="center",
                va="top",
                fontsize=6.6,
                style="italic",
                color=grey,
                zorder=4,
            )
            for j, detail in enumerate(model["detail"]):
                ax.text(
                    x + bw / 2,
                    y0 + bh - 0.185 - j * 0.052,
                    detail,
                    ha="center",
                    va="top",
                    fontsize=6.2,
                    color=grey,
                    zorder=4,
                )
            arrow(ax, x + bw / 2, y0 - 0.005, x + bw / 2, 0.505)

        ax.text(
            0.5,
            0.985,
            "Three influential models, three representational designs",
            ha="center",
            va="top",
            fontsize=7.4,
            color=grey,
        )

        # Identical geometry screen.
        box(ax, 0.0, 0.35, 1.0, 0.15, "#F2F2F2", lw=0.8)
        ax.text(
            0.5,
            0.455,
            "identical four-metric geometry screen on the gene-embedding matrix",
            ha="center",
            va="center",
            fontsize=7.2,
            fontweight="bold",
        )
        ax.text(
            0.5,
            0.395,
            "norm  $\\cdot$  centroid distance  $\\cdot$  cosine to centroid  "
            "$\\cdot$  isolation ($k=10$);   outlier if $|z|>3$ on any metric",
            ha="center",
            va="center",
            fontsize=6.4,
            color=grey,
        )

        for i, model in enumerate(models):
            x = i * (bw + gap)
            ax.text(
                x + bw / 2,
                0.305,
                model["out"],
                ha="center",
                va="center",
                fontsize=7.0,
                fontweight="bold",
            )

        # Shared tests and the model-specific downstream follow-up.
        box(ax, 0.0, 0.03, 0.665, 0.19, "#FFFFFF", lw=0.8)
        box(ax, 0.685, 0.03, 0.315, 0.19, "#FFFFFF", lw=0.8)
        ax.text(
            0.3325,
            0.175,
            "all three models",
            ha="center",
            va="center",
            fontsize=7.2,
            fontweight="bold",
        )
        ax.text(
            0.8425,
            0.175,
            "Geneformer only",
            ha="center",
            va="center",
            fontsize=7.2,
            fontweight="bold",
        )
        shared_questions = [
            "stable under\nre-calling?",
            "recur in ESM-2\nsequence space?",
            "associated with\nClinVar?",
        ]
        for i, question in enumerate(shared_questions):
            ax.text(
                0.115 + i * 0.222,
                0.085,
                question,
                ha="center",
                va="center",
                fontsize=6.3,
                color=grey,
                linespacing=1.35,
            )
        ax.text(
            0.8425,
            0.085,
            "costly to delete vs\nmatched controls?",
            ha="center",
            va="center",
            fontsize=6.3,
            color=grey,
            linespacing=1.35,
        )
        for xval in (0.2255, 0.4475):
            ax.plot([xval, xval], [0.045, 0.135], lw=0.5, color="0.75", zorder=1)

        fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)
        fig.savefig(
            REPO / "figures" / "F1_designs.pdf",
            metadata=DETERMINISTIC_PDF,
        )
    plt.close(fig)
    print("  F1_designs.pdf")


if __name__ == "__main__":
    print("Building PSB figures ->", REPO / "figures")
    fig1_designs()
    fig2_stability()
    fig3_esm2()
    fig4_nullband()
    print("Done.")
