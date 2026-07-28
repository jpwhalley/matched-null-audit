"""
E8 — Covariate-aware ClinVar association for scFM geometric outliers.

WHY THIS EXISTS
---------------
Stage 5 of the audit asks whether outlier status carries independent
disease-gene signal. An unadjusted odds ratio cannot answer that: geometric
outliers are enriched for loss-of-function-constrained genes, are highly
expressed, broadly detected, and long -- all of which independently predict
ClinVar membership. The question is whether outlier status adds anything AFTER
those covariates.

SPECIFICATION HIERARCHY -- declared here, reported in full, no cherry-picking.
This hierarchy was fixed AFTER an exploratory analysis had already been run and
three different answers seen (see the audit trail below). It is therefore NOT
pre-specified, and the manuscript must not describe it as such. The defence is
that every specification is reported, not that the primary was chosen blind.

  PRIMARY      Firth penalised logistic regression, full-coverage annotation
               table, ALL genes retained (including mitochondrial), profile-
               likelihood confidence intervals. Firth is the appropriate choice
               because all 13 MT- genes are ClinVar-positive, giving complete
               separation that breaks unpenalised maximum likelihood.
  SENSITIVITY  Ordinary logistic regression restricted to non-mitochondrial
               genes. Removes the separation by removing the separating class.
  DIAGNOSTIC   Unpenalised logistic on all genes with the mitochondrial
               covariate OMITTED. Reported ONLY to document what goes wrong
               when a perfectly separated class is left unadjusted; its
               estimates must not be quoted. This is the specification that
               produced the spurious OR 2.28.

  PRIMARY GENE SET   all Geneformer geometric outliers.
  EXPLORATORY        the Geneformer n scGPT cross-model subset. Small (n=72),
                     selected post hoc as the set with an apparent signal, and
                     therefore not a second headline test.

AUDIT TRAIL -- why this file exists at all.
An earlier ad-hoc run used the ROOT data/Table_S1.csv, which carries gene
lengths for only 11,752 of 18,915 genes; the remainder were median-imputed.
That run reported an adjusted shared-set OR of 1.95 (p=0.024). With the
full-coverage table the same model gives 1.47 (p=0.215), and retaining
mitochondrial genes gives 2.28 (p=0.007). The conclusion moved across the
significance boundary on two undocumented choices. That is precisely the
failure mode this paper is about, and it is why the analysis is scripted here
with every specification reported.

Outputs (revision/outputs/):
  E8_clinvar_adjusted.csv    -- every specification x gene set
  E8_clinvar_adjusted.json   -- primary result, universes, provenance

Usage:  python E8_clinvar_adjusted.py
"""

import hashlib
import json
import sys
import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

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
    d["clinvar"] = (d.clinvar_disease == True).astype(int)      # noqa: E712
    d["constrained"] = ((d.pLI > 0.9) | (d.LOEUF < 0.35)).astype(int)
    d["ribo"] = u.str.match(r"^(RPL|RPS|MRPL|MRPS)\d").astype(int)
    d["mito"] = u.str.startswith("MT-").astype(int)
    d["log_expr"] = np.log1p(d.max_tpm.fillna(0))
    d["breadth"] = d.expression_breadth.fillna(0)
    d["log_len"] = np.log1p(d.gene_length_bp)      # no imputation: coverage checked
    d["gf"] = d.outlier_class.isin(
        ["GF-only", "GF∩scGPT", "GF∩SF", "All three"]).astype(int)
    d["shared"] = d.outlier_class.isin(["GF∩scGPT", "All three"]).astype(int)
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
          f"({int(d.gf.sum())} Geneformer, {int(d.shared.sum())} shared)")
    print(f"  Non-mito: {len(d_nm)} "
          f"({int(d_nm.gf.sum())} Geneformer, {int(d_nm.shared.sum())} shared)")
    print()

    for tier, desc, dat, covs, fn in specs:
        print(f"  [{tier}] {desc}")
        for flag, setname in [("gf", "all Geneformer outliers"),
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

    summary = {
        "primary": {
            "specification": "Firth penalised logistic, full-coverage table, "
                             "all genes, profile-likelihood CI",
            "gene_set": "all Geneformer geometric outliers",
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
        "specification_dependence": {
            "note": "The shared-set conclusion is specification-dependent. "
                    "Report the full grid in E8_clinvar_adjusted.csv rather "
                    "than a single number.",
            "stale_table_nonmito_OR": 1.95,
            "full_table_nonmito_OR": 1.47,
            "full_table_allgene_unpenalised_OR": 2.28,
        },
        "not_pre_specified": (
            "This hierarchy was fixed after exploratory results were seen. "
            "The manuscript must not describe it as pre-specified; the "
            "defence is completeness of reporting, not blindness of choice."),
        "provenance": {
            "table": str(TABLE_S1.relative_to(REPO)),
            "sha256_16": file_hash(TABLE_S1),
            "gene_length_coverage": round(float(cov), 4),
        },
    }
    with open(OUT / "E8_clinvar_adjusted.json", "w") as f:
        json.dump(clean(summary), f, indent=2)

    print("=" * 78)
    print("  PRIMARY RESULT")
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
