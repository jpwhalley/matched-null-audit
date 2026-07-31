"""
E11 — Broad cancer-cell dependency against ClinVar, on one common universe.

EXPLORATORY. Added 2026-07-31, after the E8 result was known. Not
pre-specified, not a confirmatory test, and not to be promoted to one.

WHY THIS EXISTS
---------------
E8 finds no detectable ClinVar association for any model's outlier set. That is
a negative, and on its own it does not say what the geometry does track. This
script asks whether the outliers instead track genes that cancer cell lines
depend on.

WHAT DepMap MEASURES, AND WHAT IT DOES NOT
------------------------------------------
CRISPRGeneEffect is a Chronos gene-effect score: the effect of knocking a gene
out on proliferation and survival in cultured cancer cell lines. The scale is
anchored at 0 for reference non-essential genes and -1 for reference
common-essential genes. Averaging across lines and thresholding therefore
identifies *broad cancer-cell dependency*. It is not organismal essentiality,
not "centrality to life", and not a claim about development. The wording in
this script and in the manuscript is kept to "dependency" for that reason.

WHAT THIS SCRIPT FIXES RELATIVE TO ITS FIRST VERSION
----------------------------------------------------
1. One common universe. Genes must have complete covariates, a ClinVar call and
   a dependency score, so dependency and ClinVar are compared on identical rows.
2. Both strata. All genes, and excluding ribosomal and mitochondrial genes,
   because a class-adjusted model is inestimable (outlier ribosomal genes are
   almost all dependencies) and exclusion is the estimable alternative.
3. A continuous outcome alongside the binary ones, so nothing rests on a
   threshold.
4. Multiplicity stated explicitly rather than left implicit.
5. No claim that any association is independent of gene class. Where that model
   cannot be fitted, no number is reported.

Outputs (outputs/):
  E11_dependency.csv   -- one row per model x outcome x stratum
  E11_dependency.json  -- headline figures, caveats and provenance

Usage:
  python E11_dependency.py
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _precision import clean
from _ribosomal_panel import ribosomal_symbols, panel_provenance

_HERE = Path(__file__).resolve()
if _HERE.parent.name == "analysis":
    REPO = _HERE.parent.parent
    DATA, OUT = REPO / "data", REPO / "outputs"
else:
    REPO = _HERE.parent.parent.parent
    DATA, OUT = REPO / "notebooks" / "data", REPO / "revision" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
TABLE = DATA / "Table_S1.csv"

DEPMAP = {
    "release": "DepMap Public 25Q3",
    "file": "CRISPRGeneEffect.csv",
    "sha256": "fa44b6a7cefe8748aaa2bef21bd49aa105f7ebfe819510da3022266d88f24d51",
    "retrieved": "2026-04-14",
    "n_cell_lines": 1186,
    "n_genes": 18435,
    "portal": ("https://depmap.org/portal/data_page/?tab=allData"
               "&releasename=DepMap%20Public%2025Q3"
               "&filename=CRISPRGeneEffect.csv"),
    "statistic": "mean Chronos gene effect across cell lines",
    "measures": ("effect of knockout on proliferation and survival in cancer "
                 "cell lines; scaled to 0 for reference non-essential and -1 "
                 "for reference common-essential genes"),
    "does_not_measure": ("organismal essentiality, developmental requirement, "
                         "or importance to a whole organism"),
}
CUT, CUT_STRICT = -0.5, -1.0
COVARIATES = ["logexpr", "expression_breadth", "loglen", "constrained"]
SETS = {
    "Geneformer":   ["GF-only", "GF∩scGPT", "GF∩SF", "All three"],
    "scGPT":        ["scGPT-only", "GF∩scGPT", "All three"],
    "scFoundation": ["SF-only", "GF∩SF", "All three"],
}


def _patch_sklearn_shim():
    from sklearn.base import BaseEstimator
    if not hasattr(BaseEstimator, "_validate_data"):
        from sklearn.utils.validation import check_X_y

        def _validate_data(self, X, y=None, **kw):
            kw.pop("ensure_min_samples", None); kw.pop("dtype", None)
            return check_X_y(X, y)
        BaseEstimator._validate_data = _validate_data


def firth(d, outcome):
    _patch_sklearn_shim()
    from firthlogist import FirthLogisticRegression
    X = d[["outlier"] + COVARIATES].to_numpy(float)
    y = d[outcome].to_numpy(int)
    m = FirthLogisticRegression(max_iter=500, fit_intercept=True)
    m.fit(X, y)
    return (float(np.exp(m.coef_[0])), float(np.exp(m.ci_[0][0])),
            float(np.exp(m.ci_[0][1])), float(m.pvals_[0]), "OR")


def linear(d, outcome):
    """Continuous gene effect: coefficient is the shift in mean Chronos score.
    Negative means outliers are more strongly depended upon."""
    X = sm.add_constant(d[["outlier"] + COVARIATES].astype(float))
    m = sm.OLS(d[outcome].astype(float), X).fit()
    c, s = m.params["outlier"], m.bse["outlier"]
    return (float(c), float(c - 1.96 * s), float(c + 1.96 * s),
            float(m.pvalues["outlier"]), "beta")


def main() -> None:
    sha = hashlib.sha256(TABLE.read_bytes()).hexdigest()
    t = pd.read_csv(TABLE, low_memory=False)
    g = t.gene_symbol.astype(str).str.upper()
    t["ribo"] = g.isin(ribosomal_symbols(t)).astype(int)
    t["mito"] = g.str.startswith("MT-").astype(int)
    t["constrained"] = ((t.pLI > 0.9) | (t.LOEUF < 0.35)).astype(int)
    t["logexpr"] = np.log1p(t.max_tpm)
    t["loglen"] = np.log1p(t.gene_length_bp)
    t["clinvar"] = (t.clinvar_disease == True).astype(int)  # noqa: E712
    t["dependency"] = t.mean_dependency
    t["dep"] = (t.mean_dependency < CUT).astype(int)
    t["dep_strict"] = (t.mean_dependency < CUT_STRICT).astype(int)

    # ONE common universe: everything below is fitted on exactly these rows.
    U = t.dropna(subset=["mean_dependency", "clinvar"] + COVARIATES).copy()

    print("=" * 78)
    print("  E11 - Broad cancer-cell dependency against ClinVar  [EXPLORATORY]")
    print("=" * 78)
    print(f"  Table: {TABLE.name}  sha256[:16]={sha[:16]}")
    print(f"  Dependency: {DEPMAP['release']}, {DEPMAP['statistic']}")
    print(f"  Common universe: {len(U)} genes with complete covariates, a "
          f"ClinVar call and a dependency score")
    print(f"    dependencies (< {CUT}): {int(U.dep.sum())}   "
          f"(< {CUT_STRICT}): {int(U.dep_strict.sum())}   "
          f"ClinVar: {int(U.clinvar.sum())}")

    STRATA = [("all genes", U),
              ("excluding ribosomal and mitochondrial",
               U[(U.ribo == 0) & (U.mito == 0)])]
    OUTCOMES = [("dep", f"dependency (< {CUT})", firth),
                ("dep_strict", f"dependency (< {CUT_STRICT})", firth),
                ("dependency", "mean gene effect (continuous)", linear),
                ("clinvar", "ClinVar membership", firth)]

    rows = []
    for sname, S in STRATA:
        print(f"\n  stratum: {sname}  (n = {len(S)})")
        for ocol, olabel, fn in OUTCOMES:
            print(f"    {olabel}")
            for m, cls in SETS.items():
                d = S.copy()
                d["outlier"] = d.outlier_class.isin(cls).astype(int)
                if d.outlier.sum() < 5:
                    print(f"      {m:14s}  too few outliers in stratum")
                    continue
                try:
                    est, lo, hi, p, kind = fn(d, ocol)
                except Exception as e:
                    print(f"      {m:14s}  not estimable ({type(e).__name__})")
                    rows.append(dict(stratum=sname, outcome=ocol, model=m,
                                     estimate=None, note="not estimable"))
                    continue
                rows.append(dict(stratum=sname, outcome=ocol, model=m,
                                 n=len(d), n_outliers=int(d.outlier.sum()),
                                 estimate_kind=kind, estimate=est,
                                 ci_low=lo, ci_high=hi, p=p))
                print(f"      {m:14s}{int(d.outlier.sum()):6d}"
                      f"{est:9.3f}  ({lo:.3f}, {hi:.3f})  p={p:.1e}")

    n_tests = len(SETS) * len(OUTCOMES) * len(STRATA)
    pd.DataFrame(rows).to_csv(OUT / "E11_dependency.csv", index=False)
    (OUT / "E11_dependency.json").write_text(json.dumps({
        "status": "EXPLORATORY, not pre-specified, added 2026-07-31",
        "common_universe": int(len(U)),
        "multiplicity": (
            f"{n_tests} tests across 3 models, 4 outcomes and 2 strata. "
            f"p values are uncorrected; a Bonferroni threshold would be "
            f"{0.05 / n_tests:.4f}. Interpret accordingly."),
        "class_adjusted_model": (
            "Not fitted. Adding ribosomal and mitochondrial indicators "
            "separates, because outlier ribosomal genes are almost all "
            "dependencies. The estimable alternative is the exclusion "
            "stratum, which is reported instead. No claim is made that any "
            "association is independent of gene class."),
        "what_depmap_measures": DEPMAP["measures"],
        "what_depmap_does_not_measure": DEPMAP["does_not_measure"],
        "results": clean(rows),
        "provenance": {"table": TABLE.name, "sha256_16": sha[:16],
                       "depmap": DEPMAP, **panel_provenance()},
    }, indent=2, ensure_ascii=False))
    print(f"\n  {n_tests} tests, p values uncorrected "
          f"(Bonferroni threshold {0.05 / n_tests:.4f})")
    print(f"  Saved: {OUT / 'E11_dependency.csv'}")
    print(f"         {OUT / 'E11_dependency.json'}")


if __name__ == "__main__":
    main()
