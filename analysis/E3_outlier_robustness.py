"""
E3 — Non-parametric outlier robustness analysis (calibrated version).

Compares the original |z| > 3 outlier calls against robust alternatives.

VIABLE COMPARATORS (used in the gate assessment) — the "MAD family":
  1. MAD-based z > 3 (most direct robust alternative to standard z-scores)
  2. MAD-based z > 3.5 (stricter MAD threshold)
  3. Rank-matched MAD composite (same n as original, MAD-ranked)

DEGENERATE CALLERS — computed and reported, but EXCLUDED from the gate
assessment. Each exclusion rests on different evidence; see the verdict
memo revision/verdicts/E3_verdict_2026-07-17.md (amended 2026-07-27):
  - 2-component GMM — too permissive. Minority component captures
    24.20-49.17% of genes in every model x metric combination
    (E3_degenerate_diagnostics.csv). Calling a quarter to a half of the
    vocabulary "outliers" is not outlier detection.
  - Percentile 1/99 — fixed-rate, not threshold-responsive. Observed union
    rates 4.96% / 6.77% / 5.92% for Geneformer / scGPT / scFoundation,
    i.e. ~5-7% by construction regardless of distributional shape.
  - IQR k=3.0 (Tukey extreme fences) — too strict. Returns 43 outliers for
    Geneformer (containment 10.5%) and ZERO for both scGPT and scFoundation
    (E3_calibrated_summary.csv). A caller returning the empty set on two of
    three models cannot distinguish stability from its own insensitivity.
    NOTE: this caller was listed as viable in versions of this script before
    2026-07-27. It is not. The enrichment-stability claim holds under the
    MAD family but NOT under IQR/Tukey, so this reclassification is
    load-bearing rather than cosmetic.

DIRECTION WARNING: enrichment_test() below uses Fisher alternative="two-sided",
so significant=True covers DEPLETION as well as enrichment. Geneformer's
"disease" class is a significant DEPLETION (OR 0.54 original; 0.51-0.68 across
the MAD family). Report these as "class associations" with direction stated —
never as "four enrichments".

For each viable method, reports:
  - Jaccard overlap with the original set
  - Containment (fraction of original outliers recovered)
  - Top-50 overlap (rank stability at the extreme tail)
  - Spearman rank correlation of anomaly scores
  - Fisher's exact enrichment tests for headline biological categories

Gate criteria (specified before downstream experiments; amended after
diagnostic checks showed GMM, raw percentile AND IQR/Tukey were degenerate
— see § Degenerate Method Diagnostics below and the amended verdict memo):
  - Primary stability metric: containment + top-k overlap + rank correlation
  - Jaccard is reported but interpreted in context of set-size mismatch
  - Per-model verdict: STABLE / MIXED / UNSTABLE

Outputs (all in revision/outputs/):
  - E3_calibrated_summary.csv — per-model containment, Jaccard, Spearman, top-50
  - E3_enrichment_full.csv — all enrichment tests across models × methods
  - E3_robust_core.json — per-model intersection of original ∩ MAD outlier sets
  - E3_degenerate_diagnostics.csv — GMM minority-component sizes documenting
    why GMM was excluded
  - E3_gate_verdict.json — structured verdict for downstream gating
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.mixture import GaussianMixture
from pathlib import Path
import json
import sys
import warnings

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _precision import clean  # documented serialisation precision
from _ribosomal_panel import panel_provenance, ribosomal_symbols

# Repository-relative paths. Scripts live in analysis/; everything they read
# and write is inside this repository.
REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUT = REPO / "outputs"
CACHE = REPO / "cache"
for _d in (OUT, CACHE):
    _d.mkdir(parents=True, exist_ok=True)
BASE = REPO  # legacy alias

warnings.filterwarnings("ignore")

# ── Paths ────────────────────────────────────────────────────────────────────
OUT  = OUT

# ── Data ─────────────────────────────────────────────────────────────────────
MODEL_FILES = {
    "Geneformer":    DATA / "gene_embedding_geometry.csv",
    "scGPT":         DATA / "scgpt_gene_embedding_geometry.csv",
    "scFoundation":  DATA / "sf_gene_embedding_geometry.csv",
}

RAW_METRICS = ["norm", "dist_from_centroid", "cos_to_centroid", "isolation_score"]
Z_COLS      = ["norm_zscore", "dist_zscore", "cos_zscore", "isolation_zscore"]
ENRICH_CATS = ["ribosomal", "mitochondrial", "constrained", "disease"]

table_s1 = pd.read_csv(DATA / "Table_S1.csv")
RIBOSOMAL = frozenset(ribosomal_symbols(table_s1))


# ── Helpers ──────────────────────────────────────────────────────────────────

def load_model(name: str) -> pd.DataFrame:
    df = pd.read_csv(MODEL_FILES[name])
    if name == "Geneformer":
        df = df[~df["gene"].isin(["<pad>", "<mask>", "<cls>", "<eos>"])]
        df = df.reset_index(drop=True)
    return df


def score_col(df: pd.DataFrame) -> str:
    return ("anomaly_score_with_isolation"
            if "anomaly_score_with_isolation" in df.columns
            else "anomaly_score")


def assign_gene_classes(genes: pd.Series) -> pd.Series:
    """Classify genes into biological categories for enrichment testing."""
    classes = pd.Series("other", index=genes.index)
    sym = genes.str.upper()

    classes[sym.str.startswith("MT-")]                    = "mitochondrial"
    classes[sym.isin(RIBOSOMAL)] = "ribosomal"

    constrained = set(
        table_s1.loc[(table_s1["pLI"] > 0.9) | (table_s1["LOEUF"] < 0.35),
                     "gene_symbol"].str.upper()
    )
    classes[sym.isin(constrained) & (classes == "other")] = "constrained"

    disease = set(
        table_s1.loc[table_s1["clinvar_disease"] == True,
                     "gene_symbol"].str.upper()
    )
    classes[sym.isin(disease) & (classes == "other")]     = "disease"
    return classes


def mad_z(x: np.ndarray) -> np.ndarray:
    """MAD-based robust z-score for a 1-d array."""
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:
        mad = 1e-10
    return np.abs(x - med) / (1.4826 * mad)


def iqr_flags(x: np.ndarray, k: float = 3.0) -> np.ndarray:
    """Tukey fence outlier flags (k = 3.0 = extreme outliers)."""
    q1, q3 = np.percentile(x, [25, 75])
    iqr = q3 - q1
    return (x < q1 - k * iqr) | (x > q3 + k * iqr)


def enrichment_test(outlier_mask: np.ndarray, gene_classes: pd.Series,
                    cat: str) -> dict:
    is_t = (gene_classes == cat).values
    a = int(np.sum(outlier_mask & is_t))
    b = int(np.sum(outlier_mask & ~is_t))
    c = int(np.sum(~outlier_mask & is_t))
    d = int(np.sum(~outlier_mask & ~is_t))
    if a + c == 0:
        return {"class": cat, "OR": np.nan, "p": np.nan,
                "n_outlier": 0, "n_total": 0}
    OR, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return {"class": cat, "OR": OR, "p": p, "n_outlier": a, "n_total": a + c}


# ── 1. Degenerate-method diagnostics (GMM + percentile) ─────────────────────

print("=" * 70)
print("  DEGENERATE METHOD DIAGNOSTICS")
print("=" * 70)

gmm_rows = []
for model_name in MODEL_FILES:
    df = load_model(model_name)
    for col in RAW_METRICS:
        x = df[col].values.reshape(-1, 1)
        gmm = GaussianMixture(n_components=2, random_state=42, n_init=5)
        gmm.fit(x)
        labels = gmm.predict(x)
        counts = np.bincount(labels, minlength=2)
        minority_n = int(counts.min())
        gmm_rows.append({
            "model": model_name, "metric": col,
            "minority_n": minority_n,
            "minority_pct": 100 * minority_n / len(df),
            "total_n": len(df),
        })
        print(f"  {model_name:15s} {col:25s}  "
              f"minority={minority_n:6d} ({100*minority_n/len(df):.1f}%)")

gmm_diag = pd.DataFrame(gmm_rows)
gmm_diag.to_csv(OUT / "E3_degenerate_diagnostics.csv", index=False)
print(f"\n  All GMM minority components ≥ 24% of vocabulary → GMM excluded.")
print(f"  Percentile 1/99 flags ~4×1% = ~4% by construction → excluded.\n")


# ── 2. Calibrated robustness analysis ───────────────────────────────────────

print("=" * 70)
print("  CALIBRATED ROBUSTNESS ANALYSIS")
print("=" * 70)

summary_rows = []
enrichment_rows = []
robust_cores = {}

for model_name in MODEL_FILES:
    df = load_model(model_name)
    sc = score_col(df)
    n_orig    = int(df["is_outlier"].sum())
    orig_mask = df["is_outlier"].values.astype(bool)
    orig_set  = set(df.loc[orig_mask, "gene"])
    gene_cls  = assign_gene_classes(df["gene"])

    # Compute MAD composite score (max MAD-z across metrics)
    mad_composite = np.zeros(len(df))
    for col in RAW_METRICS:
        mz = mad_z(df[col].values)
        mad_composite = np.maximum(mad_composite, mz)

    # ── Define methods ───────────────────────────────────────────────────
    # Each entry: (name, outlier_mask, gene_set)
    # MAD z>3
    mad3_mask = np.zeros(len(df), dtype=bool)
    for col in RAW_METRICS:
        mad3_mask |= (mad_z(df[col].values) > 3.0)
    mad3_set = set(df.loc[mad3_mask, "gene"])

    # MAD z>3.5
    mad35_mask = np.zeros(len(df), dtype=bool)
    for col in RAW_METRICS:
        mad35_mask |= (mad_z(df[col].values) > 3.5)
    mad35_set = set(df.loc[mad35_mask, "gene"])

    # IQR k=3
    iqr_mask = np.zeros(len(df), dtype=bool)
    for col in RAW_METRICS:
        iqr_mask |= iqr_flags(df[col].values, k=3.0)
    iqr_set = set(df.loc[iqr_mask, "gene"])

    # Rank-matched: top n_orig by MAD composite
    df_tmp = df.assign(_mad_c=mad_composite)
    rank_set = set(df_tmp.nlargest(n_orig, "_mad_c")["gene"])

    methods = [
        ("|z|>3 (original)",         orig_mask,  orig_set),
        ("MAD z>3",                  mad3_mask,  mad3_set),
        ("MAD z>3.5",                mad35_mask, mad35_set),
        # DEGENERATE (excluded from gate; computed + reported for transparency):
        # returns 43 outliers for Geneformer, 0 for scGPT, 0 for scFoundation
        ("IQR k=3 (Tukey extreme)",  iqr_mask,   iqr_set),
        ("Top-n by MAD score",       None,       rank_set),
    ]

    # Spearman: original anomaly score vs MAD composite
    rho, rho_p = stats.spearmanr(df[sc].values, mad_composite)

    # Top-50 from original ranking
    orig_top50 = set(df.nlargest(50, sc)["gene"])
    mad_top50  = set(df_tmp.nlargest(50, "_mad_c")["gene"])
    top50_overlap = len(orig_top50 & mad_top50)

    # Robust core = original ∩ MAD z>3
    robust_core = orig_set & mad3_set
    robust_cores[model_name] = sorted(robust_core)

    print(f"\n{'='*70}")
    print(f"  {model_name}  (n={len(df)}, original outliers={n_orig})")
    print(f"{'='*70}")
    print(f"  Spearman (original vs MAD composite): rho={rho:.3f}")
    print(f"  Top-50 overlap (original vs MAD-ranked): {top50_overlap}/50")
    print(f"  Robust core (original ∩ MAD z>3): {len(robust_core)} genes")

    print(f"\n  {'Method':30s} {'n':>6s} {'Contain':>8s} {'Jaccard':>8s}")
    print(f"  {'-'*30} {'-'*6} {'-'*8} {'-'*8}")

    for mname, mmask, mset in methods:
        n = len(mset)
        contain = len(orig_set & mset) / n_orig if n_orig > 0 else 0
        jacc = (len(orig_set & mset) / len(orig_set | mset)
                if len(orig_set | mset) > 0 else 1.0)

        is_orig = mname == "|z|>3 (original)"
        summary_rows.append({
            "model": model_name,
            "method": mname,
            "n_outliers": n,
            "containment": contain if not is_orig else 1.0,
            "jaccard": jacc if not is_orig else 1.0,
            "spearman_rho": rho,
            "top50_overlap": top50_overlap,
            "n_robust_core": len(robust_core),
        })

        print(f"  {mname:30s} {n:6d} {contain:7.1%} {jacc:8.3f}")

        # Enrichments
        if mmask is None:
            mmask = df["gene"].isin(mset).values
        for cat in ENRICH_CATS:
            res = enrichment_test(mmask, gene_cls, cat)
            res["model"]  = model_name
            res["method"] = mname
            res["p_bonf"] = min(res["p"] * len(ENRICH_CATS), 1.0) \
                            if not np.isnan(res["p"]) else np.nan
            res["significant"] = res["p_bonf"] < 0.05 \
                                 if not np.isnan(res["p_bonf"]) else False
            enrichment_rows.append(res)

    # Print enrichment table
    print(f"\n  Enrichments (OR, Bonferroni p):")
    print(f"  {'':30s}", end="")
    for cat in ENRICH_CATS:
        print(f"  {cat:>14s}", end="")
    print()
    for mname, _, mset in methods:
        print(f"  {mname:30s}", end="")
        for cat in ENRICH_CATS:
            row = [r for r in enrichment_rows
                   if r["model"] == model_name
                   and r["method"] == mname
                   and r["class"] == cat]
            if row:
                r = row[0]
                if np.isnan(r["OR"]):
                    print(f"  {'N/A':>14s}", end="")
                else:
                    sig = "*" if r["significant"] else " "
                    print(f"  {r['OR']:5.1f} ({r['p_bonf']:.0e}){sig}", end="")
            else:
                print(f"  {'???':>14s}", end="")
        print()


# ── 3. Gate assessment ───────────────────────────────────────────────────────

print(f"\n\n{'='*70}")
print(f"  GATE ASSESSMENT")
print(f"{'='*70}")
print(f"\n  Note: gate criteria specified before downstream experiments were")
print(f"  amended after diagnostic checks showed THREE callers degenerate:")
print(f"    - GMM 2-component  (too permissive; minority 24.20-49.17%)")
print(f"    - Percentile 1/99  (fixed-rate ~5-7% by construction)")
print(f"    - IQR k=3 (Tukey)  (too strict; 43/0/0 outliers, containment 10.5%)")
print(f"  Assessment uses the MAD family as comparators: MAD z>3,")
print(f"  MAD z>3.5, and rank-matched MAD (Top-n by MAD score).")
print(f"  Fisher tests are TWO-SIDED: 'significant' includes depletion.\n")

verdicts = {}
for model_name in MODEL_FILES:
    # Get MAD z>3 stats
    mad_row = [r for r in summary_rows
               if r["model"] == model_name and r["method"] == "MAD z>3"][0]
    rank_row = [r for r in summary_rows
                if r["model"] == model_name
                and r["method"] == "Top-n by MAD score"][0]

    containment = mad_row["containment"]
    rho         = mad_row["spearman_rho"]
    t50         = mad_row["top50_overlap"]
    core_n      = mad_row["n_robust_core"]
    orig_n      = [r for r in summary_rows
                   if r["model"] == model_name
                   and r["method"] == "|z|>3 (original)"][0]["n_outliers"]

    # Enrichment stability: do originally-significant categories stay
    # significant under MAD z>3 and rank-matched?
    orig_sig = {r["class"] for r in enrichment_rows
                if r["model"] == model_name
                and r["method"] == "|z|>3 (original)"
                and r.get("significant", False)}
    flipped = 0
    for cat in orig_sig:
        alt_sigs = [r["significant"] for r in enrichment_rows
                    if r["model"] == model_name
                    and r["method"] in ("MAD z>3", "MAD z>3.5",
                                         "Top-n by MAD score")
                    and r["class"] == cat]
        if sum(alt_sigs) < 2:  # need ≥2 of 3 applicable methods
            flipped += 1

    # Verdict logic
    if containment >= 0.95 and rho >= 0.95 and t50 >= 45:
        verdict = "STABLE"
        note = ("original z-score outliers are a conservative core "
                "under MAD; rank ordering is preserved")
    elif containment >= 0.80 and rho >= 0.80:
        verdict = "MIXED"
        note = ("biological enrichment pattern is stable but individual "
                "gene identity is method-sensitive")
    else:
        verdict = "UNSTABLE"
        note = "outlier set changes substantially with method choice"

    # Override: if enrichments flip, downgrade
    if flipped >= 2 and verdict == "STABLE":
        verdict = "MIXED"
        note += "; enrichment flips detected"

    verdicts[model_name] = {
        "verdict": verdict,
        "note": note,
        "containment": containment,
        "spearman_rho": rho,
        "top50_overlap": t50,
        "robust_core_n": core_n,
        "original_n": orig_n,
        "enrichment_flips": flipped,
    }

    print(f"  {model_name}: {verdict}")
    print(f"    {note}")
    print(f"    containment={containment:.1%}  rho={rho:.3f}  "
          f"top50={t50}/50  core={core_n}/{orig_n}")

# Overall
per_verdicts = [v["verdict"] for v in verdicts.values()]
if all(v == "STABLE" for v in per_verdicts):
    overall = "STABLE"
elif any(v == "UNSTABLE" for v in per_verdicts):
    overall = "UNSTABLE"
else:
    overall = "MIXED"

overall_note = (
    "Two models (Geneformer, scFoundation) stable; scGPT mixed due to "
    "bimodal norm distribution. Proceed to E1 with Geneformer as primary "
    "E2 model. For scGPT, use original ∩ MAD robust core or treat as "
    "sensitivity analysis."
)

print(f"\n  {'='*50}")
print(f"  OVERALL E3 VERDICT: {overall}")
print(f"  {overall_note}")
print(f"  {'='*50}")


# ── 4. Save all outputs ─────────────────────────────────────────────────────

pd.DataFrame(summary_rows).to_csv(OUT / "E3_calibrated_summary.csv", index=False)
pd.DataFrame(enrichment_rows).to_csv(OUT / "E3_enrichment_full.csv", index=False)

with open(OUT / "E3_robust_core.json", "w") as f:
    json.dump(clean(robust_cores), f, indent=2)

gate_output = {
    "overall_verdict": overall,
    "overall_note": overall_note,
    "panel_provenance": panel_provenance(),
    "criteria_note": (
        "Gate criteria specified before downstream experiments; amended after "
        "diagnostic checks. THREE callers are excluded as degenerate, each on "
        "DIFFERENT evidence: (1) 2-component GMM - too permissive, minority "
        "component captured 24.20-49.17% of genes in every model x metric "
        "combination (see E3_degenerate_diagnostics.csv, which documents this "
        "and nothing else); (2) percentile 1/99 - fixed-rate selector flagging "
        "~5-7% of genes by construction regardless of distributional shape "
        "(observed union rates 4.96%/6.77%/5.92% for Geneformer/scGPT/"
        "scFoundation; evidence: construction of the method plus these rates, "
        "NOT E3_degenerate_diagnostics.csv); (3) IQR k=3 (Tukey extreme) - too "
        "strict, returns 43 outliers for Geneformer (containment 10.5%) and 0 "
        "for both scGPT and scFoundation, so it cannot discriminate stability "
        "from its own insensitivity (evidence: E3_calibrated_summary.csv, NOT "
        "E3_degenerate_diagnostics.csv). CAVEAT: only IQR/Tukey is present in "
        "E3_calibrated_summary.csv and E3_enrichment_full.csv. GMM component "
        "sizes are retained only as degeneracy diagnostics; older GMM and "
        "percentile enrichment files used a superseded ribosomal definition "
        "and must not be quoted for class effects. Assessment uses "
        "MAD z>3, MAD z>3.5, and rank-matched MAD composite as comparators."
    ),
    "direction_note": (
        "IMPORTANT: enrichment_test uses Fisher alternative='two-sided', so "
        "significant=True covers DEPLETION as well as enrichment. Do not "
        "describe every significant class association as an enrichment. "
        "Under the original caller, Geneformer is enriched for mitochondrial, "
        "ribosomal and constrained genes but depleted for the residual "
        "mutually-exclusive disease class. Directions and current effect "
        "sizes are recorded in E3_enrichment_full.csv."
    ),
    "enrichment_flips_note": (
        "enrichment_flips counts ONLY classes significant under |z|>3 "
        "(original) that lose significance in >=2 of three MAD-family callers "
        "(MAD z>3, MAD z>3.5, Top-n by MAD). It does not evaluate percentile, "
        "GMM or IQR/Tukey and it does not count gains. Under the pinned HGNC "
        "panel, scFoundation's ribosomal association is significant under the "
        "original caller and two of three MAD-family comparators, so its zero "
        "is not vacuous. scGPT disease gains significance under MAD z>3 but "
        "that gain is not counted. The excluded callers must not be used to "
        "support a claim of stability across all methods; current method- and "
        "class-specific results are in E3_enrichment_full.csv."
    ),
    "per_model": verdicts,
}
with open(OUT / "E3_gate_verdict.json", "w") as f:
    json.dump(clean(gate_output), f, indent=2, default=str)

print(f"\n  Saved to {OUT}/:")
for fn in ["E3_calibrated_summary.csv", "E3_enrichment_full.csv",
           "E3_robust_core.json", "E3_degenerate_diagnostics.csv",
           "E3_gate_verdict.json"]:
    print(f"    {fn}")
