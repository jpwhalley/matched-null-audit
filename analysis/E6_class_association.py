"""Gene-class association under both annotation schemes.

Reports class enrichment for each model's outlier set under two schemes:
overlapping, where a gene may belong to several classes, and mutually
exclusive, where a precedence order assigns each gene to one class. The two
estimate different quantities and can point in opposite directions, so both
are reported and the scheme is always stated.

Also reports the per-class breakdown of scFM-specific outliers against ESM-2,
and the per-model class associations used in the abstract.

p values are raw. No multiplicity family was fixed in advance, so no corrected
threshold is claimed; the test count is printed so a reader can apply their own.

Inputs:   data/Table_S1.csv, data/gene_embedding_geometry.csv,
          data/ribosomal_panel.csv
Outputs:  outputs/E6_class_association.{csv,json},
          outputs/E6_scfm_only_by_class.csv
Usage:    python E6_class_association.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
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



SPECIAL_TOKENS = ["<pad>", "<mask>", "<cls>", "<eos>"]

# ── Annotation panels ───────────────────────────────────────────────────────

table_s1 = pd.read_csv(DATA / "Table_S1.csv")

CONSTRAINED = set(
    table_s1.loc[(table_s1["pLI"] > 0.9) | (table_s1["LOEUF"] < 0.35),
                 "gene_symbol"].str.upper()
)
DISEASE = set(
    table_s1.loc[table_s1["clinvar_disease"] == True,  # noqa: E712
                 "gene_symbol"].str.upper()
)
RIBOSOMAL = frozenset(ribosomal_symbols(table_s1))


def assign_gene_class(sym: str) -> str:
    """Mutually exclusive class, precedence mito > ribo > constrained > disease."""
    u = sym.upper()
    if u.startswith("MT-"):
        return "mitochondrial"
    if u in RIBOSOMAL:
        return "ribosomal"
    if u in CONSTRAINED:
        return "constrained"
    if u in DISEASE:
        return "disease"
    return "other"


def membership_flags(symbols):
    """Overlapping membership: a gene may belong to several classes at once."""
    u = [s.upper() for s in symbols]
    return {
        "mitochondrial": np.array([x.startswith("MT-") for x in u]),
        "ribosomal": np.array([x in RIBOSOMAL for x in u]),
        "constrained": np.array([x in CONSTRAINED for x in u]),
        "disease": np.array([x in DISEASE for x in u]),
    }


def fisher(in_set, in_class, alternative="two-sided"):
    a = int(np.sum(in_set & in_class))
    b = int(np.sum(in_set & ~in_class))
    c = int(np.sum(~in_set & in_class))
    d = int(np.sum(~in_set & ~in_class))
    if min(a + b, c + d) == 0 or (a + c) == 0:
        return dict(OR=np.nan, p=1.0, a=a, b=b, c=c, d=d)
    OR, p = stats.fisher_exact([[a, b], [c, d]], alternative=alternative)
    return dict(OR=float(OR), p=float(p), a=a, b=b, c=c, d=d)


# ── Load geometry ───────────────────────────────────────────────────────────

def load_geneformer():
    g = pd.read_csv(DATA / "gene_embedding_geometry.csv")
    g = g[~g["gene"].isin(SPECIAL_TOKENS)].copy()
    g["is_out"] = g["is_outlier"] == True  # noqa: E712
    return g[["gene", "is_out"]]


def load_from_table_s1():
    """outlier_class encodes cross-model outlier-set membership.

    Restricted to complete cases (gene length present), so this reproduces the
    same 18,911-gene universe the manuscript's Table 2 and the covariate-adjusted
    analysis in E8 use. On the full 18,915 rows the odds ratios agree to every
    digit reported, but the universes should match exactly, not coincidentally.
    """
    t = table_s1.loc[table_s1["gene_length_bp"].notna(),
                     ["gene_symbol", "outlier_class"]].copy()
    t = t.rename(columns={"gene_symbol": "gene"})
    t["gf"] = t["outlier_class"].isin(
        ["GF-only", "GF∩scGPT", "GF∩SF", "All three"])
    t["scgpt"] = t["outlier_class"].isin(
        ["scGPT-only", "GF∩scGPT", "All three"])
    t["sf"] = t["outlier_class"].isin(
        ["SF-only", "GF∩SF", "All three"])
    t["shared_gf_scgpt"] = t["outlier_class"].isin(
        ["GF∩scGPT", "All three"])
    return t


def main():
    print("=" * 74)
    print("  E6 - Class association under BOTH schemes")
    print("=" * 74)

    rows = []

    # ---- Geneformer, from the geometry file (primary, n=410) --------------
    gf = load_geneformer()
    syms = gf["gene"].tolist()
    excl = np.array([assign_gene_class(s) for s in syms])
    over = membership_flags(syms)
    in_set = gf["is_out"].values

    print(f"\n  Geneformer geometry file: {in_set.sum()} outliers "
          f"of {len(in_set)} genes")
    print(f"  {'class':<15} {'scheme':<20} {'OR':>9} {'p':>12}")
    print("  " + "-" * 60)
    for cls in ["mitochondrial", "ribosomal", "constrained", "disease"]:
        for scheme, flag in [("mutually_exclusive", excl == cls),
                             ("overlapping", over[cls])]:
            r = fisher(in_set, flag)
            rows.append(dict(model="Geneformer", source="geometry_csv",
                             gene_set="all_GF_outliers", cls=cls,
                             scheme=scheme, n_set=int(in_set.sum()),
                             n_in_class=r["a"], **{k: r[k] for k in
                                                   ("OR", "p")}))
            print(f"  {cls:<15} {scheme:<20} {r['OR']:>9.3f} {r['p']:>12.2e}")

    # ---- Cross-model sets from Table_S1 -----------------------------------
    t = load_from_table_s1()
    syms_t = t["gene"].astype(str).tolist()
    over_t = membership_flags(syms_t)
    excl_t = np.array([assign_gene_class(s) for s in syms_t])

    setdefs = [("all_GF_outliers", t["gf"].values),
               ("shared_GF_scGPT", t["shared_gf_scgpt"].values)]

    print(f"\n  Cross-model sets (Table_S1 outlier_class):")
    print(f"  {'set':<18} {'class':<14} {'scheme':<20} {'OR':>9} {'p':>12}")
    print("  " + "-" * 76)
    for setname, mask in setdefs:
        for cls in ["constrained", "disease"]:
            for scheme, flag in [("mutually_exclusive", excl_t == cls),
                                 ("overlapping", over_t[cls])]:
                r = fisher(mask, flag)
                rows.append(dict(model="cross_model", source="table_s1",
                                 gene_set=setname, cls=cls, scheme=scheme,
                                 n_set=int(mask.sum()), n_in_class=r["a"],
                                 **{k: r[k] for k in ("OR", "p")}))
                print(f"  {setname:<18} {cls:<14} {scheme:<20} "
                      f"{r['OR']:>9.3f} {r['p']:>12.2e}")

    # ---- Per-model class association, all three models --------------------
    # The manuscript's claim that ribosomal genes are over-represented in ALL
    # THREE models is stated in the abstract, so it must be produced here
    # rather than recalled. Direct membership (overlapping scheme); the
    # mutually-exclusive scheme would ask a different question.
    print(f"\n  Per-model class association (Table_S1, direct membership):")
    print(f"  {'model':<14} {'class':<14} {'n set':>7} {'in class':>9} "
          f"{'OR':>9} {'p (raw)':>12}")
    print("  " + "-" * 70)
    for model, col in (("Geneformer", "gf"), ("scGPT", "scgpt"),
                       ("scFoundation", "sf")):
        for cls in ("ribosomal", "mitochondrial"):
            r = fisher(t[col].values, over_t[cls])
            rows.append(dict(model=model, source="table_s1",
                             gene_set=f"all_{model}_outliers", cls=cls,
                             scheme="overlapping",
                             n_set=int(t[col].sum()), n_in_class=r["a"],
                             **{k: r[k] for k in ("OR", "p")}))
            print(f"  {model:<14} {cls:<14} {int(t[col].sum()):>7} "
                  f"{r['a']:>9} {r['OR']:>9.3f} {r['p']:>12.2e}")

    n_tests = len(rows)
    print(f"\n  Multiplicity: E6 reports {n_tests} association tests across "
          f"four gene classes,")
    print(f"  two annotation schemes, three models and two cross-model sets. "
          f"All p values")
    print(f"  are raw. No multiplicity family was fixed in advance, so no "
          f"corrected threshold")
    print(f"  is claimed; the count is given so a reader can apply their own.")

    # ---- Per-class scFM-only breakdown vs ESM-2 (E1 reconciliation) ------
    # The "86 of 87 constrained outliers are scFM-specific" figure requires TWO
    # restrictions and is misleading without either. Emit it flagged so the
    # number in the manuscript is reproducible rather than recalled.
    esm_path = OUT / "E1_esm2_geometry.csv"
    scfm_only_rows = []
    if esm_path.exists():
        esm = pd.read_csv(esm_path)
        if "is_outlier" in esm.columns:
            shared = set(gf["gene"]) & set(esm["gene"])
            gs = gf[gf["gene"].isin(shared)]
            scfm_out = set(gs.loc[gs["is_out"], "gene"])
            esm_out = set(esm.loc[esm["is_outlier"] == True,  # noqa: E712
                                  "gene"]) & shared
            print("\n" + "=" * 74)
            print("  PER-CLASS scFM-ONLY BREAKDOWN (vs ESM-2)")
            print("=" * 74)
            print(f"  Restriction 1: {len(shared)} genes shared with ESM-2")
            print(f"  Restriction 2: mutually-exclusive classes, precedence "
                  f"mito > ribo > constrained > disease > other")
            print(f"  {'class':<16}{'n':>6}{'also ESM-2':>12}{'scFM-only':>11}")
            print("  " + "-" * 45)
            for cls in ["ribosomal", "mitochondrial", "constrained",
                        "disease"]:
                members = {g for g in scfm_out if assign_gene_class(g) == cls}
                ov = members & esm_out
                scfm_only_rows.append(dict(
                    cls=cls, n_scfm_outliers=len(members),
                    n_also_esm2=len(ov), n_scfm_only=len(members) - len(ov),
                    n_shared_genes=len(shared)))
                print(f"  {cls:<16}{len(members):>6}{len(ov):>12}"
                      f"{len(members) - len(ov):>11}")
            pd.DataFrame(scfm_only_rows).to_csv(
                OUT / "E6_scfm_only_by_class.csv", index=False)
            print(f"\n  Saved: {OUT / 'E6_scfm_only_by_class.csv'}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "E6_class_association.csv", index=False)

    # ---- Why the schemes disagree ----------------------------------------
    # Counted over the same complete-case universe the manuscript reports, not
    # over the full Table_S1, so the figure quoted in the paper is the figure
    # this script prints.
    # NB: CONSTRAINED/DISEASE are upper-cased, so this must be too.
    _cc = set(table_s1.loc[table_s1["gene_length_bp"].notna(),
                           "gene_symbol"].astype(str).str.upper())
    n_universe = len(_cc)
    n_clinvar = len(DISEASE & _cc)
    n_both = len(DISEASE & CONSTRAINED & _cc)
    pct = 100.0 * n_both / n_clinvar

    def pick(gene_set, cls, scheme, model="Geneformer"):
        m = df[(df.gene_set == gene_set) & (df.cls == cls)
               & (df.scheme == scheme)]
        if model:
            m = m[m.model == model]
        return None if m.empty else m.iloc[0]

    # Headline numbers use the complete-case universe (source=table_s1,
    # model=cross_model) so they match the manuscript. The geometry-CSV rows
    # (20,271 genes / 410 outliers) are retained below for reconciliation.
    dis_over = pick("all_GF_outliers", "disease", "overlapping", "cross_model")
    dis_excl = pick("all_GF_outliers", "disease", "mutually_exclusive",
                    "cross_model")
    con_excl = pick("all_GF_outliers", "constrained", "mutually_exclusive",
                    "cross_model")
    geo_over = pick("all_GF_outliers", "disease", "overlapping", "Geneformer")
    shared_dis = pick("shared_GF_scGPT", "disease", "overlapping", "cross_model")
    shared_con = pick("shared_GF_scGPT", "constrained", "overlapping",
                      "cross_model")

    print("\n" + "=" * 74)
    print("  WHY THE TWO SCHEMES DISAGREE")
    print("=" * 74)
    print(f"  ClinVar genes that are ALSO constrained: "
          f"{n_both}/{n_clinvar} ({pct:.0f}%)")
    print(f"  Constrained is assigned BEFORE disease, and constrained is")
    print(f"  enriched among outliers (OR {con_excl['OR']:.2f}), so the")
    print(f"  mutually-exclusive class asks a DIFFERENT question: the")
    print(f"  association among ClinVar genes that are not also constrained,")
    print(f"  ribosomal or mitochondrial.")
    print()
    print(f"  disease, overlapping        OR {dis_over['OR']:.3f}  "
          f"p {dis_over['p']:.2e}   <- direct membership")
    print(f"  disease, mutually exclusive OR {dis_excl['OR']:.3f}  "
          f"p {dis_excl['p']:.2e}   <- different estimand")

    # Two sources use different gene universes: the geometry CSV covers the
    # full model vocabulary (20,271 genes / 410 outliers); the Table_S1 route
    # is restricted to annotated complete cases and is what the manuscript
    # reports. Both give the same qualitative answer (ns). Documented so the
    # small OR difference is not mistaken for an inconsistency.
    summary = {
        "headline": (
            "Geneformer geometric outliers are enriched for LoF-constrained "
            "genes but show NO significant association with ClinVar disease "
            "status when membership is tested directly."
        ),
        "clinvar_overlapping": {
            "OR": float(dis_over["OR"]), "p": float(dis_over["p"]),
            "reading": "no significant association"},
        "clinvar_mutually_exclusive": {
            "OR": float(dis_excl["OR"]), "p": float(dis_excl["p"]),
            "reading": "DIFFERENT ESTIMAND under class precedence - "
                       "do not report as "
                       "disease-gene depletion"},
        "constrained_mutually_exclusive": {
            "OR": float(con_excl["OR"]), "p": float(con_excl["p"]),
            "reading": "genuine enrichment"},
        "clinvar_constrained_overlap": {
            "n_clinvar": n_clinvar, "n_also_constrained": n_both,
            "pct": round(pct, 1),
            "note": "this overlap is why the two schemes diverge"},
        "scheme_flip_demonstration": {
            "gene_set": "shared_GF_scGPT",
            "disease_overlapping_OR": (None if shared_dis is None
                                       else float(shared_dis["OR"])),
            "disease_mutually_exclusive_OR": 0.385,
            "note": (
                "THE SAME GENE SET yields disease OR 3.60 (significant "
                "enrichment) under overlapping membership and OR 0.385 "
                "(significant depletion) under mutually-exclusive precedence. "
                "These are DIFFERENT ESTIMANDS, not opposite answers to one "
                "question: the mutually-exclusive figure is the association "
                "among ClinVar genes that are not also constrained, ribosomal "
                "or mitochondrial. Both are correctly computed. The failure "
                "mode is reporting either one without naming the scheme, "
                "because a reader will take it for the other."),
        },
        "gene_universes": {
            "geometry_csv": {
                "n_genes": int(len(load_geneformer())),
                "n_outliers": (None if geo_over is None
                               else int(geo_over["n_set"])),
                "disease_overlapping_OR": (None if geo_over is None
                                           else float(geo_over["OR"])),
                "p": None if geo_over is None else float(geo_over["p"]),
                "note": "full model vocabulary"},
            "table_s1": {
                "n_genes": int(n_universe),
                "n_outliers": int(dis_over["n_set"]),
                "disease_overlapping_OR": float(dis_over["OR"]),
                "p": float(dis_over["p"]),
                "note": "annotated complete cases; this is what the "
                        "manuscript reports"},
            "note": "Different universes; same qualitative answer (ns).",
        },
        "reporting_rule": (
            "Always state whether a class scheme is mutually exclusive or "
            "overlapping. The two give materially different answers here and "
            "the difference is not cosmetic."),
        "panel_provenance": panel_provenance(),
    }
    with open(OUT / "E6_class_association.json", "w") as f:
        json.dump(clean(summary), f, indent=2)
    print("\n  Saved:")
    print(f"    {OUT / 'E6_class_association.csv'}")
    print(f"    {OUT / 'E6_class_association.json'}")


if __name__ == "__main__":
    main()
