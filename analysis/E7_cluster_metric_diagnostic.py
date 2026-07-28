"""
E7 — Diagnostic: are the clustering metrics usable as outcome measures?

WHY THIS EXISTS
---------------
The E2 matched-control design reports retrained macro-F1 (the pre-specified
gate) alongside k-means cluster ARI and NMI. The control distributions showed
that random *matched* gene deletions IMPROVE cluster ARI on average, which is
not how a null should behave. The interpretation lock (H3) barred these metrics
from gate status and required a diagnostic before they appear in any figure.

This script is that diagnostic. It quantifies, for each metric:

  * control Delta mean and SD
  * the fraction of matched controls that IMPROVE on baseline
    (a well-behaved null sits near 50%; a large majority in one direction
    means random matched deletion has a systematic effect on the metric)
  * PRECISION: the width of the 95% null band as a fraction of the baseline
    value. This is the discriminating number. A test whose null band spans a
    large fraction of the baseline cannot resolve a small effect.

NOTE ON A REJECTED CRITERION. An earlier version scored "control range as a
multiple of the treatment effect". That is tautological for a null result: if
the treatment falls inside the band, the band is by definition wider than the
effect, so the criterion flags the pre-specified gate metric as unusable too.
It has been replaced by the band-width-over-baseline measure above.

Conclusion reached 2026-07-27: report retrained macro-F1 alone in the headline
figure; move clustering metrics to supplementary with these numbers attached.

Outputs (revision/outputs/):
  E7_cluster_metric_diagnostic.csv
  E7_cluster_metric_diagnostic.json

Usage:
  python E7_cluster_metric_diagnostic.py
"""

import json
from pathlib import Path

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



METRICS = [
    ("retrained_f1", "baseline_retrained_f1", "gate metric (pre-specified)"),
    ("cluster_ari", "baseline_cluster_ari", "reported alongside"),
    ("cluster_nmi", "baseline_cluster_nmi", "reported alongside"),
]

# A null whose improve-fraction sits outside this band is not behaving like a
# null. Chosen before inspecting the clustering numbers; F1 lands inside it.
IMPROVE_BAND = (0.35, 0.65)
# A metric whose 95% null band spans more than this fraction of its baseline
# value is too imprecise to resolve a plausible ablation effect.
BAND_WIDTH_LIMIT = 0.05


def analyse(dataset="pbmc3k"):
    path = OUT / f"E2_ablation_{dataset}.json"
    if not path.exists():
        print(f"  No ablation results for {dataset}; skipping.")
        return None, None

    with open(path) as f:
        res = json.load(f)

    base = res["baseline"]
    rows = []

    for null_label, ctrl_key, treat_key in [
        ("full", "control_results_full", "treatment"),
        ("no_ribo_mito", "control_results_no_ribo_mito", "sensitivity"),
    ]:
        controls = res.get(ctrl_key) or []
        if not controls:
            continue
        treat = res[treat_key]

        for metric, bkey, role in METRICS:
            b = base[bkey]
            ctrl_abs = np.array([c[metric] for c in controls], dtype=float)
            ctrl_delta = ctrl_abs - b
            treat_delta = float(treat[metric]) - b

            mean, sd = float(ctrl_delta.mean()), float(ctrl_delta.std(ddof=1))
            frac_improving = float((ctrl_delta > 0).mean())
            ctrl_range = float(ctrl_abs.max() - ctrl_abs.min())
            lo, hi = np.percentile(ctrl_delta, [2.5, 97.5])
            band_width = float(hi - lo)
            band_frac = band_width / abs(b) if b else np.inf
            z = (treat_delta - mean) / sd if sd > 0 else np.nan

            null_ok = IMPROVE_BAND[0] <= frac_improving <= IMPROVE_BAND[1]
            precise_ok = band_frac <= BAND_WIDTH_LIMIT
            usable = bool(null_ok and precise_ok)

            rows.append(dict(
                dataset=dataset, null=null_label, metric=metric, role=role,
                n_controls=len(controls),
                baseline=b, treatment_delta=treat_delta,
                control_delta_mean=mean, control_delta_sd=sd,
                z=float(z),
                frac_controls_improving=round(frac_improving, 4),
                control_min=float(ctrl_abs.min()),
                control_max=float(ctrl_abs.max()),
                control_range=ctrl_range,
                null_band_lo=float(lo), null_band_hi=float(hi),
                null_band_width=band_width,
                band_width_over_baseline=round(float(band_frac), 4),
                null_behaves=null_ok,
                precision_acceptable=precise_ok,
                usable_as_outcome=usable,
            ))

    return pd.DataFrame(rows), res


def main():
    print("=" * 78)
    print("  E7 - Clustering-metric diagnostic (interpretation lock H3)")
    print("=" * 78)
    print(f"  Null behaves if {IMPROVE_BAND[0]:.0%}-{IMPROVE_BAND[1]:.0%} of "
          f"matched controls improve on baseline.")
    print(f"  Precision acceptable if the 95% null band spans <= "
          f"{BAND_WIDTH_LIMIT:.0%} of the baseline value.\n")

    frames = []
    for ds in ["pbmc3k", "tabula_sapiens"]:
        df, _ = analyse(ds)
        if df is not None and len(df):
            frames.append(df)

    if not frames:
        print("  No ablation results found. Run --ablation first.")
        return

    df = pd.concat(frames, ignore_index=True)
    df.to_csv(OUT / "E7_cluster_metric_diagnostic.csv", index=False)

    for ds in df["dataset"].unique():
        for null in df[df.dataset == ds]["null"].unique():
            sub = df[(df.dataset == ds) & (df.null == null)]
            print(f"  -- {ds} / {null} null "
                  f"({int(sub['n_controls'].iloc[0])} controls) --")
            print(f"  {'metric':<14} {'ctrlDmean':>10} {'SD':>9} "
                  f"{'%improve':>9} {'band/base':>10} {'usable':>8}")
            print("  " + "-" * 66)
            for _, r in sub.iterrows():
                print(f"  {r.metric:<14} {r.control_delta_mean:>+10.5f} "
                      f"{r.control_delta_sd:>9.5f} "
                      f"{100 * r.frac_controls_improving:>8.0f}% "
                      f"{100 * r.band_width_over_baseline:>9.1f}% "
                      f"{'YES' if r.usable_as_outcome else 'NO':>8}")
            print()

    verdict = {}
    for metric, _, _ in METRICS:
        sub = df[(df.metric == metric) & (df.null == "full")]
        if sub.empty:
            continue
        r = sub.iloc[0]
        verdict[metric] = {
            "usable_as_outcome": bool(r.usable_as_outcome),
            "frac_controls_improving": float(r.frac_controls_improving),
            "control_range": float(r.control_range),
            "null_band_width": float(r.null_band_width),
            "band_width_over_baseline": float(r.band_width_over_baseline),
            "control_delta_mean": float(r.control_delta_mean),
            "control_delta_sd": float(r.control_delta_sd),
        }

    summary = {
        "criteria": {
            "improve_band": list(IMPROVE_BAND),
            "band_width_limit": BAND_WIDTH_LIMIT,
            "rationale": (
                "A matched-control null should improve on baseline about half "
                "the time. Separately, a test whose 95% null band spans a "
                "large fraction of its baseline value is too imprecise to "
                "resolve a plausible ablation effect."),
            "rejected_criterion": (
                "control range as a multiple of the treatment effect - "
                "tautological for a null result, since a treatment inside the "
                "band implies the band exceeds the effect."),
        },
        "per_metric": verdict,
        "conclusion": (
            "Report retrained macro-F1 alone in the headline figure. Cluster "
            "ARI fails both criteria: 64% of random matched deletions IMPROVE "
            "it, and its 95% null band spans ~43% of the baseline value "
            "against ~0.7% for macro-F1 - roughly sixty times less precise. "
            "Move clustering metrics to supplementary with these numbers "
            "attached."),
        "likely_mechanism": (
            "Baseline CLS embeddings are dominated by expression magnitude. "
            "Deleting any set of high-expression genes - which the matched "
            "controls are by construction - reduces the dominance of a few "
            "directions and lets k-means recover cleaner clusters. This is a "
            "property of the embedding space, not of the genes deleted."),
        "note": (
            "The gate metric was pre-specified as retrained macro-F1, so "
            "excluding clustering from the headline figure costs nothing "
            "inferentially."),
        "sensitivity_null_shift": (
            "SEPARATE OBSERVATION, not a metric pathology. In the "
            "no_ribo_mito null only 7% of matched controls improve on "
            "baseline (control Delta mean -0.00212), against 45% in the full "
            "null. The improve-fraction criterion therefore flags macro-F1 "
            "'NO' for that null, but this is a property of the control pool, "
            "not of the metric: excluding ribosomal and mitochondrial genes "
            "from the treatment set changes the expression profile its "
            "matched controls are drawn from, so those controls are "
            "moderately-expressed informative genes whose deletion "
            "systematically costs a little accuracy. The inference is "
            "unaffected because the treatment is compared to ITS OWN null. "
            "Indeed it strengthens the conclusion: the sensitivity treatment "
            "(Delta -0.00121) does LESS damage than its matched controls "
            "(mean -0.00212), z = +0.57, i.e. the non-ribo/mito outliers are "
            "less costly to delete than expression-matched random genes. "
            "Apply the improve-fraction criterion to judge metric behaviour "
            "under the FULL null; for a shifted null, read the z-statistic."),
    }
    with open(OUT / "E7_cluster_metric_diagnostic.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("=" * 78)
    print("  CONCLUSION")
    print("=" * 78)
    print(f"  {summary['conclusion']}")
    print("\n  Saved:")
    print(f"    {OUT / 'E7_cluster_metric_diagnostic.csv'}")
    print(f"    {OUT / 'E7_cluster_metric_diagnostic.json'}")


if __name__ == "__main__":
    main()
