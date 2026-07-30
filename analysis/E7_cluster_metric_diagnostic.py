"""
E7 — Diagnostic: are the clustering metrics usable as outcome measures?

WHY THIS EXISTS
---------------
The E2 matched-control design reports retrained macro-F1 alongside k-means
cluster ARI and NMI. The dated plan listed annotation accuracy and clustering
stability as outcomes but did not designate a gate metric; macro-F1 was chosen
after that plan was written. The interpretation lock (H3) barred the clustering
metrics from gate status unless a diagnostic showed that their control
distributions were precise enough to resolve an effect of the size at issue.

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
effect, so the criterion flags the retained gate metric as unusable too.
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
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _precision import clean  # documented serialisation precision

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
    ("retrained_f1", "baseline_retrained_f1", "headline gate"),
    ("cluster_ari", "baseline_cluster_ari", "pre-specified secondary outcome"),
    ("cluster_nmi", "baseline_cluster_nmi", "pre-specified secondary outcome"),
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

    # Prefer the pinned standalone baseline, as E2 evaluate() and Figure 4 do.
    # The probe baseline moves ~6e-4 between environments and the ablation
    # stores whichever value was in memory when it ran.
    #
    # This matters for the DESCRIPTIVE improve-fraction, which moves by more
    # than ten points between the two baselines -- precisely why that quantity
    # is no longer a usability gate. The
    # surviving criterion, band width over baseline, is insensitive to a shift
    # this small: the band width is unchanged and the denominator moves by
    # 0.06%.
    base = dict(res["baseline"])
    baseline_source = "ablation-run"
    pinned_path = OUT / f"E2_baseline_{dataset}.json"
    if pinned_path.exists():
        with open(pinned_path) as f:
            pin = json.load(f)
        for k in ("baseline_retrained_f1", "baseline_cluster_ari",
                  "baseline_cluster_nmi"):
            if k in pin:
                base[k] = pin[k]
        baseline_source = pinned_path.name
        if not pin.get("environment"):
            print(f"  WARNING: {pinned_path.name} carries no environment "
                  f"fingerprint.")
    print(f"  Baseline source: {baseline_source}")
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

            # improve-fraction is DESCRIPTIVE ONLY -- see the header note.
            null_ok = IMPROVE_BAND[0] <= frac_improving <= IMPROVE_BAND[1]
            precise_ok = band_frac <= BAND_WIDTH_LIMIT
            usable = bool(precise_ok)

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
    print(f"  USABILITY: the 95% null band must span <= "
          f"{BAND_WIDTH_LIMIT:.0%} of the baseline value.")
    print(f"  %improve is reported for description only and does not affect "
          f"the verdict.\n")

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
            "usability_rests_on": "band_width_over_baseline only",
            "improve_fraction_status": (
                "DESCRIPTIVE ONLY. Not a usability gate. It is sensitive to "
                "the baseline at a scale far below the metric's own, so at "
                "some baselines it would reject the retained gate metric. "
                "Live values are in improve_fraction_observed and the CSV."),
            "improve_fraction_observed": {
                f"{r.dataset}/{r.null}/{r.metric}":
                    round(float(r.frac_controls_improving), 4)
                for r in df.itertuples()},
            "rationale": (
                "A test whose 95% null band spans a large fraction of its "
                "baseline value is too imprecise to resolve a plausible "
                "ablation effect. That is the sole usability criterion."),
            "rejected_criterion": (
                "control range as a multiple of the treatment effect - "
                "tautological for a null result, since a treatment inside the "
                "band implies the band exceeds the effect."),
        },
        "per_metric": verdict,
        "conclusion": (
            "Report retrained macro-F1 alone in the headline figure. Cluster "
            "ARI is excluded on PRECISION: its 95% null band spans ~42% of the "
            "baseline value against ~0.7% for macro-F1, roughly sixty times "
            "less precise, so it cannot resolve an effect of the size at "
            "issue. The earlier 'random deletion reliably improves ARI' "
            "argument does NOT survive the corrected controls and is not "
            "relied on. Move clustering metrics to supplementary with these "
            "numbers attached."),
        "possible_mechanism_not_tested": (
            "The broad clustering nulls could reflect sensitivity of k-means "
            "to dominant embedding directions, but this mechanism was not "
            "tested and no interpretation depends on it."),
        "note": (
            "The dated plan listed annotation accuracy and clustering stability "
            "but did not designate a gate metric. Macro-F1 was selected after "
            "the plan and before the Tabula Sapiens result was analysed. The "
            "clustering exclusion is reported rather than treated as "
            "pre-specified, and no headline inference rests on it."),
        "sensitivity_null_shift": (
            "SEPARATE OBSERVATION, not a metric pathology. The no_ribo_mito "
            "control pool differs from the full pool: excluding ribosomal and "
            "mitochondrial genes from the treatment set changes the expression "
            "profile its matched controls are drawn from, so those controls "
            "are moderately-expressed informative genes whose deletion "
            "systematically costs a little accuracy. The inference is "
            "unaffected because the treatment is compared to ITS OWN null; "
            "read the z-statistic, not the improve-fraction. Numbers for both "
            "nulls are in the CSV and are regenerated with the controls."),
    }
    with open(OUT / "E7_cluster_metric_diagnostic.json", "w") as f:
        json.dump(clean(summary), f, indent=2)

    print("=" * 78)
    print("  CONCLUSION")
    print("=" * 78)
    print(f"  {summary['conclusion']}")
    print("\n  Saved:")
    print(f"    {OUT / 'E7_cluster_metric_diagnostic.csv'}")
    print(f"    {OUT / 'E7_cluster_metric_diagnostic.json'}")


if __name__ == "__main__":
    main()
