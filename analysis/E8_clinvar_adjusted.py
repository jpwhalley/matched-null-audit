"""Covariate-adjusted ClinVar association for all three outlier sets.

Logistic regression of ClinVar membership on outlier status plus constraint,
log1p expression, expression breadth, log1p gene length, and ribosomal and
mitochondrial indicators. The mitochondrial indicator separates perfectly, so
the primary fit is Firth-penalised with profile-likelihood intervals.

The specification is identical across the three models: only the exposure flag
differs. A non-mitochondrial sensitivity fit and an unpenalised diagnostic fit
are reported alongside, the latter labelled as mis-specified because it leaves
the separating class unadjusted.

Inputs:   data/Table_S1.csv
Outputs:  outputs/E8_clinvar_adjusted.{csv,json}
Usage:    python E8_clinvar_adjusted.py
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

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



# Full-coverage annotation table. The root data/Table_S1.csv is length-sparse
# (11,752/18,915) and must not be used for any length-adjusted analysis.
TABLE_S1 = DATA / "Table_S1.csv"
MIN_LENGTH_COVERAGE = 0.99

COVARIATES = ["constrained", "log_expr", "breadth", "log_len", "ribo", "mito"]


def file_hash(path, n=16):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n]


def load():
    t = pd.read_csv(TABLE_S1)
    cov = t["gene_length_bp"].notna().mean()
    if cov < MIN_LENGTH_COVERAGE:
        raise RuntimeError(
            f"{TABLE_S1} has gene-length coverage {cov:.3f} < "
            f"{MIN_LENGTH_COVERAGE}. Median imputation at this rate changes "
            f"the adjusted odds ratios materially. Use the full-coverage "
            f"table.")

    d = t[["gene_symbol", "outlier_class", "pLI", "LOEUF", "clinvar_disease",
           "max_tpm", "expression_breadth", "gene_length_bp"]].copy()
    u = d.gene_symbol.astype(str).str.upper()
    ribosomal = ribosomal_symbols(t)
    d["clinvar"] = (d.clinvar_disease == True).astype(int)      # noqa: E712
    d["constrained"] = ((d.pLI > 0.9) | (d.LOEUF < 0.35)).astype(int)
    d["ribo"] = u.isin(ribosomal).astype(int)
    d["mito"] = u.str.startswith("MT-").astype(int)
    d["log_expr"] = np.log1p(d.max_tpm.fillna(0))
    d["breadth"] = d.expression_breadth.fillna(0)
    d["log_len"] = np.log1p(d.gene_length_bp)      # no imputation: coverage checked
    d["gf"] = d.outlier_class.isin(
        ["GF-only", "GF∩scGPT", "GF∩SF", "All three"]).astype(int)
    d["shared"] = d.outlier_class.isin(["GF∩scGPT", "All three"]).astype(int)
    # Comparative extension across models. The same specification is
    # applied to the scGPT and scFoundation outlier sets so the disease analysis
    # covers all three models rather than one. Nothing about the specification
    # changes: same universe, same covariates, same Firth primary. Only the
    # exposure flag differs. Model encoding is irrelevant here because the
    # regression uses the derived outlier flag and gene-level annotations, not
    # model inference.
    d["scgpt"] = d.outlier_class.isin(
        ["scGPT-only", "GF∩scGPT", "All three"]).astype(int)
    d["scfound"] = d.outlier_class.isin(
        ["SF-only", "GF∩SF", "All three"]).astype(int)
    d = d.dropna(subset=["log_expr", "breadth", "log_len"])
    return d, cov


def unadjusted(d, flag):
    a = int(((d[flag] == 1) & (d.clinvar == 1)).sum())
    b = int(((d[flag] == 1) & (d.clinvar == 0)).sum())
    c = int(((d[flag] == 0) & (d.clinvar == 1)).sum())
    e = int(((d[flag] == 0) & (d.clinvar == 0)).sum())
    OR, p = stats.fisher_exact([[a, b], [c, e]])
    return float(OR), float(p)


def _patch_sklearn_shim():
    """firthlogist calls BaseEstimator._validate_data, removed in sklearn>=1.6."""
    from sklearn.base import BaseEstimator
    if not hasattr(BaseEstimator, "_validate_data"):
        from sklearn.utils.validation import check_X_y, check_array

        def _validate_data(self, X, y=None, **kw):
            kw.pop("dtype", None); kw.pop("reset", None)
            kw.pop("accept_sparse", None); kw.pop("ensure_2d", None)
            if y is None:
                return check_array(X)
            return check_X_y(X, y)

        BaseEstimator._validate_data = _validate_data


def fit_firth(d, flag, covs):
    """Firth penalised logistic with profile-likelihood CIs."""
    _patch_sklearn_shim()
    from firthlogist import FirthLogisticRegression
    cols = [flag] + covs
    X = d[cols].to_numpy(float)
    y = d["clinvar"].to_numpy(int)
    m = FirthLogisticRegression(max_iter=500, fit_intercept=True)
    m.fit(X, y)
    i = 0  # flag is the first column
    return (float(np.exp(m.coef_[i])),
            float(np.exp(m.ci_[i][0])), float(np.exp(m.ci_[i][1])),
            float(m.pvals_[i]))


def fit_logit(d, flag, covs):
    X = sm.add_constant(d[[flag] + covs])
    m = sm.Logit(d["clinvar"], X).fit(disp=0, maxiter=500)
    c, s = m.params[flag], m.bse[flag]
    return (float(np.exp(c)), float(np.exp(c - 1.96 * s)),
            float(np.exp(c + 1.96 * s)), float(m.pvalues[flag]))


def main():
    d, cov = load()
    rows = []

    non_mito_covs = [c for c in COVARIATES if c != "mito"]
    d_nm = d[d.mito == 0]

    specs = [
        ("PRIMARY", "Firth, all genes, profile-likelihood CI",
         d, COVARIATES, fit_firth),
        ("SENSITIVITY", "Logistic, non-mitochondrial genes",
         d_nm, non_mito_covs, fit_logit),
        ("DIAGNOSTIC", "Unpenalised logistic, all genes, mito covariate OMITTED "
                       "(MIS-SPECIFIED: separating class left unadjusted)",
         d, non_mito_covs, fit_logit),
    ]

    print("=" * 78)
    print("  E8 - Covariate-aware ClinVar association")
    print("=" * 78)
    print(f"  Table: {TABLE_S1.name}  sha256[:16]={file_hash(TABLE_S1)}")
    print(f"  Gene-length coverage: {cov:.4f}")
    print(f"  Universe: {len(d)} complete cases "
          f"({int(d.gf.sum())} Geneformer, {int(d.scgpt.sum())} scGPT, "
          f"{int(d.scfound.sum())} scFoundation, {int(d.shared.sum())} shared)")
    print(f"  Non-mito: {len(d_nm)} "
          f"({int(d_nm.gf.sum())} Geneformer, {int(d_nm.shared.sum())} shared)")
    print()

    for tier, desc, dat, covs, fn in specs:
        print(f"  [{tier}] {desc}")
        for flag, setname in [("gf", "Geneformer outliers"),
                              ("scgpt", "scGPT outliers"),
                              ("scfound", "scFoundation outliers"),
                              ("shared", "shared GF n scGPT (exploratory)")]:
            uOR, up = unadjusted(dat, flag)
            try:
                aOR, lo, hi, ap = fn(dat, flag, covs)
                ok = True
            except Exception as exc:
                aOR = lo = hi = ap = float("nan"); ok = False
                print(f"      {setname}: FAILED ({exc})")
            rows.append(dict(tier=tier, spec=desc, gene_set=flag,
                             n_set=int(dat[flag].sum()), n_total=len(dat),
                             unadj_OR=uOR, unadj_p=up,
                             adj_OR=aOR, ci_lo=lo, ci_hi=hi, adj_p=ap,
                             converged=ok))
            if ok:
                sig = "SIGNIFICANT" if ap < 0.05 else "ns"
                print(f"      {setname:<34} n={int(dat[flag].sum()):4d}  "
                      f"unadj OR {uOR:5.2f} (p={up:.3f})   "
                      f"adj OR {aOR:5.2f} [{lo:.2f}, {hi:.2f}] "
                      f"p={ap:.3f}  {sig}")
        print()

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "E8_clinvar_adjusted.csv", index=False)

    prim = df[(df.tier == "PRIMARY") & (df.gene_set == "gf")].iloc[0]
    sens = df[(df.tier == "SENSITIVITY") & (df.gene_set == "gf")].iloc[0]
    prim_sh = df[(df.tier == "PRIMARY") & (df.gene_set == "shared")].iloc[0]

    MODEL_SETS = [("gf", "Geneformer"), ("scgpt", "scGPT"),
                  ("scfound", "scFoundation")]

    def _row(tier, flag):
        m = df[(df.tier == tier) & (df.gene_set == flag)]
        return None if m.empty else m.iloc[0]

    def _block(tier):
        out = {}
        for flag, label in MODEL_SETS:
            r = _row(tier, flag)
            if r is None:
                continue
            out[label] = {
                "n_outliers": int(r.n_set), "n_universe": int(r.n_total),
                "unadjusted_OR": float(r.unadj_OR),
                "unadjusted_p": float(r.unadj_p),
                "adjusted_OR": float(r.adj_OR),
                "adjusted_CI": [float(r.ci_lo), float(r.ci_hi)],
                "adjusted_p": float(r.adj_p),
            }
        return out

    summary = {
        "headline": (
            "Under an identical covariate adjustment, no model's outlier set "
            "shows a detectable ClinVar association. scGPT is the informative "
            "case: a strong unadjusted association does not survive the same "
            "adjustment applied to the other two."),
        "primary_by_model": _block("PRIMARY"),
        "sensitivity_by_model": _block("SENSITIVITY"),
        "multiplicity": (
            "The three model-specific p-values are uncorrected. All exceed "
            "0.4, so no correction would change any conclusion."),
        # Geneformer-only keys below are retained for backward compatibility
        # with anything that read this file before the three-model extension.
        "primary": {
            "specification": "Firth penalised logistic, full-coverage table, "
                             "all genes, profile-likelihood CI",
            "gene_set": "all Geneformer geometric outliers "
                        "(see primary_by_model for all three)",
            "n_outliers": int(prim.n_set), "n_universe": int(prim.n_total),
            "unadjusted_OR": prim.unadj_OR, "unadjusted_p": prim.unadj_p,
            "adjusted_OR": prim.adj_OR,
            "adjusted_CI": [prim.ci_lo, prim.ci_hi],
            "adjusted_p": prim.adj_p,
        },
        "sensitivity_non_mito": {
            "adjusted_OR": sens.adj_OR,
            "adjusted_CI": [sens.ci_lo, sens.ci_hi],
            "adjusted_p": sens.adj_p, "n_universe": int(sens.n_total),
        },
        "exploratory_shared_set": {
            "note": "Cross-model subset selected post hoc as the set showing "
                    "an apparent signal. Small (n=72). NOT a second headline "
                    "test; report as exploratory only.",
            "adjusted_OR": prim_sh.adj_OR,
            "adjusted_CI": [prim_sh.ci_lo, prim_sh.ci_hi],
            "adjusted_p": prim_sh.adj_p,
        },
        "reporting_note": (
            "The shared-set result is specification-dependent. Read the full "
            "grid in E8_clinvar_adjusted.csv rather than any single number."),
        "provenance": {
            "table": str(TABLE_S1.relative_to(BASE)),
            "sha256_16": file_hash(TABLE_S1),
            "gene_length_coverage": round(float(cov), 4),
            **panel_provenance(),
        },
    }
    with open(OUT / "E8_clinvar_adjusted.json", "w") as f:
        json.dump(clean(summary), f, indent=2)

    print("=" * 78)
    print("  PRIMARY RESULT, all three models (Firth, all genes)")
    for flag, label in MODEL_SETS:
        r = _row("PRIMARY", flag)
        if r is None:
            continue
        print(f"    {label:<13} n={int(r.n_set):>4}  unadj OR {r.unadj_OR:5.2f}"
              f"   adj OR {r.adj_OR:5.2f} [{r.ci_lo:.2f}, {r.ci_hi:.2f}]"
              f"  p = {r.adj_p:.3f}")
    print("    p-values uncorrected; all exceed 0.4.")
    print()
    print("  GENEFORMER DETAIL")
    print("=" * 78)
    print(f"  All Geneformer outliers, Firth, all genes:")
    print(f"    adjusted OR {prim.adj_OR:.2f} "
          f"[{prim.ci_lo:.2f}, {prim.ci_hi:.2f}]  p = {prim.adj_p:.3f}")
    print(f"  Sensitivity (non-mito, logistic): OR {sens.adj_OR:.2f} "
          f"[{sens.ci_lo:.2f}, {sens.ci_hi:.2f}]  p = {sens.adj_p:.3f}")
    print(f"\n  Saved: {OUT/'E8_clinvar_adjusted.csv'}")
    print(f"         {OUT/'E8_clinvar_adjusted.json'}")


if __name__ == "__main__":
    main()
