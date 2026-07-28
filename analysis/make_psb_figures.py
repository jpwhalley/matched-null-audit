"""
Build PSB 2027 figures F2-F4 from saved analysis outputs.

Design constraints (12-page limit, four figures):
  * each figure must fit ~0.45-0.55 page => single-column width 4.6in,
    height 2.0-2.6in
  * greyscale-safe: distinguish by shape/hatch/position, not colour alone
  * no chartjunk; the null band in F4 must be unambiguous at a glance

F1 (audit schematic + metric distributions) is built separately from the
original figure pipeline.

Outputs: psb2027/figs/{F2_stability,F3_esm2,F4_nullband}.pdf
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
def fig4_nullband(dataset="pbmc3k"):
    """THE figure. Treatment delta vs matched-control null. macro-F1 only."""
    path = OUT / f"E2_ablation_{dataset}.json"
    if not path.exists():
        print(f"  (skipping F4: no {path.name})")
        return
    res = json.load(open(path))

    canon_p = OUT / f"E2_baseline_{dataset}.json"
    base = (json.load(open(canon_p))["baseline_retrained_f1"]
            if canon_p.exists()
            else res["baseline"]["baseline_retrained_f1"])

    panels = [("Full treatment\n(50 genes)", "control_results_full",
               "treatment"),
              ("Sensitivity\n(36 genes, no ribo/mito)",
               "control_results_no_ribo_mito", "sensitivity")]

    fig, axes = plt.subplots(1, 2, figsize=(W, 2.15), sharex=True)
    for ax, (title, ck, tk) in zip(axes, panels):
        ctrls = res.get(ck) or []
        if not ctrls:
            ax.set_visible(False)
            continue
        cd = np.array([c["retrained_f1"] for c in ctrls]) - base
        td = res[tk]["retrained_f1"] - base
        lo, hi = np.percentile(cd, [2.5, 97.5])
        z = (td - cd.mean()) / cd.std(ddof=1)

        ax.axvspan(lo, hi, color="0.88", zorder=0)
        ax.hist(cd, bins=26, color=GREY, edgecolor="black", linewidth=0.35,
                zorder=2)
        ax.axvline(td, color=ACC, lw=1.6, zorder=4)
        ax.axvline(lo, color=DARK, lw=0.8, ls="--", zorder=3)
        ax.axvline(hi, color=DARK, lw=0.8, ls="--", zorder=3)

        ymax = ax.get_ylim()[1]
        ax.annotate("treatment", xy=(td, ymax * 0.94),
                    xytext=(td - abs(td) * 0.9 - 0.0016, ymax * 0.94),
                    fontsize=7, ha="right", va="center",
                    arrowprops=dict(arrowstyle="->", lw=0.7))
        ax.text(0.03, 0.86, f"$z$ = {z:+.2f}\nINSIDE", transform=ax.transAxes,
                fontsize=7, va="top",
                bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="0.6",
                          lw=0.5))
        ax.set_title(title, fontsize=7.5)
        ax.set_xlabel("$\\Delta$ macro-$F_1$ vs baseline")
    axes[0].set_ylabel(f"matched control sets\n(n = {len(res['control_results_full'])})")
    fig.tight_layout(pad=0.35)
    fig.savefig(REPO / "figures" / f"F4_nullband_{dataset}.pdf", metadata=DETERMINISTIC_PDF)
    plt.close(fig)
    print(f"  F4_nullband_{dataset}.pdf")


# ---------------------------------------------------------------- F1
def fig1_workflow():
    """The five-stage claim-auditing workflow, with this paper's outcomes.

    For a workflow paper the schematic earns the space more than metric
    distributions do: it shows that each stage can falsify independently, and
    that our claim survives two stages and fails two.
    """
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    stages = [
        ("1", "Geometry\nscreen", "which genes are\nextreme?", "410 outliers", "n"),
        ("2", "Caller\nrobustness", "artefact of\nthresholding?", "survives", "p"),
        ("3", "ESM-2\ncontrol", "generic\nbiology?", "model-specific", "p"),
        ("4", "Matched\ndeletion", "does it affect\na task?", "inside null", "f"),
        ("5", "Annotation\nscheme", "robust to\nconvention?", "reverses", "f"),
    ]

    fig, ax = plt.subplots(figsize=(W, 1.72))
    ax.set_xlim(0, 10); ax.set_ylim(0.15, 3.0); ax.axis("off")

    bw, gap = 1.72, 0.30
    x0 = 0.12
    for i, (num, name, question, outcome, kind) in enumerate(stages):
        x = x0 + i * (bw + gap)
        face = {"n": "0.93", "p": "0.86", "f": "0.62"}[kind]
        ax.add_patch(FancyBboxPatch(
            (x, 1.30), bw, 0.92, boxstyle="round,pad=0.045",
            facecolor=face, edgecolor="black", linewidth=0.7))
        ax.text(x + bw / 2, 2.02, num, ha="center", va="center",
                fontsize=6.5, color="0.35")
        ax.text(x + bw / 2, 1.70, name, ha="center", va="center",
                fontsize=7.2, linespacing=1.15)
        ax.text(x + bw / 2, 1.06, question, ha="center", va="top",
                fontsize=6.1, color="0.32", linespacing=1.2, style="italic")
        # outcome chip
        mark = {"n": "", "p": "✓ ", "f": "✗ "}[kind]
        ax.text(x + bw / 2, 2.52, mark + outcome, ha="center", va="center",
                fontsize=6.6,
                bbox=dict(boxstyle="round,pad=0.22",
                          fc="white" if kind != "f" else "0.90",
                          ec="0.45", lw=0.55))
        if i < len(stages) - 1:
            ax.add_patch(FancyArrowPatch(
                (x + bw + 0.03, 1.76), (x + bw + gap - 0.03, 1.76),
                arrowstyle="-|>", mutation_scale=7, lw=0.7, color="0.3"))

    # legend: swatch then label, laid out horizontally so nothing overlaps
    ax.add_patch(FancyBboxPatch((0.12, 0.30), 0.26, 0.17,
                                boxstyle="round,pad=0.02", facecolor="0.86",
                                edgecolor="black", linewidth=0.5))
    ax.text(0.50, 0.385, "claim survives this stage", fontsize=6.4,
            color="0.3", va="center")
    ax.add_patch(FancyBboxPatch((3.15, 0.30), 0.26, 0.17,
                                boxstyle="round,pad=0.02", facecolor="0.62",
                                edgecolor="black", linewidth=0.5))
    ax.text(3.53, 0.385, "claim fails this stage", fontsize=6.4,
            color="0.3", va="center")
    fig.tight_layout(pad=0.2)
    fig.savefig(REPO / "figures" / "F1_workflow.pdf", metadata=DETERMINISTIC_PDF)
    plt.close(fig)
    print("  F1_workflow.pdf")


if __name__ == "__main__":
    print("Building PSB figures ->", REPO / "figures")
    fig1_workflow()
    fig2_stability()
    fig3_esm2()
    fig4_nullband("pbmc3k")
    fig4_nullband("tabula_sapiens")
    print("Done.")
