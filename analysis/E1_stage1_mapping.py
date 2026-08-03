"""
E1 Stage 1 — Gene-to-protein mapping feasibility for ESM-2 comparison.

Maps each model's gene vocabulary to a canonical protein using mygene.info
(backed by NCBI/Ensembl/UniProt). Produces a feasibility table and bias
audit to determine whether the ESM-2 comparison (Stage 2) can proceed
without structural bias.

Gate criteria:
  PROCEED = ≥75% protein-coding outliers map, Fisher OR < 3.0
  SCOPE-RESTRICT = 50-75% mapping or OR 3.0-5.0
  STOP = <50% mapping or OR > 5.0

Outputs (all in outputs/):
  E1_feasibility_table.csv — one row per (model, gene)
  E1_bias_audit.csv — per-model mapping rates and bias statistics
  E1_mapping_cache.json — raw mygene results for reproducibility
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import json
import mygene
import sys
import warnings
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
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

OUT  = OUT

# Load Table_S1 for gene class assignment
table_s1 = pd.read_csv(DATA / "Table_S1.csv")
RIBOSOMAL = frozenset(ribosomal_symbols(table_s1))
CONSTRAINED = frozenset(
    table_s1.loc[(table_s1["pLI"] > 0.9) | (table_s1["LOEUF"] < 0.35),
                 "gene_symbol"].astype(str).str.upper()
)
DISEASE = frozenset(
    table_s1.loc[table_s1["clinvar_disease"] == True,  # noqa: E712
                 "gene_symbol"].astype(str).str.upper()
)

# Load E3 robust core for outlier definitions
with open(OUT / "E3_robust_core.json") as f:
    robust_cores = json.load(f)


# ── Gene class assignment (same as E3) ───────────────────────────────────────

def assign_gene_class(sym: str) -> str:
    sym_u = str(sym).upper()
    if sym_u.startswith("MT-"):
        return "mitochondrial"
    if sym_u in RIBOSOMAL:
        return "ribosomal"
    if sym_u in CONSTRAINED:
        return "constrained"
    if sym_u in DISEASE:
        return "disease"
    return "other"


# ── Load model vocabularies ──────────────────────────────────────────────────

print("Loading model vocabularies...")

# Geneformer: has Ensembl IDs
gf = pd.read_csv(DATA / "gene_embedding_geometry.csv")
gf = gf[~gf["gene"].isin(["<pad>", "<mask>", "<cls>", "<eos>"])].reset_index(drop=True)

# scGPT: gene symbols only
sg = pd.read_csv(DATA / "scgpt_gene_embedding_geometry.csv")

# scFoundation: gene symbols
sf = pd.read_csv(DATA / "sf_gene_embedding_geometry.csv")

print(f"  Geneformer: {len(gf)} genes (Ensembl IDs)")
print(f"  scGPT:      {len(sg)} genes (symbols)")
print(f"  scFoundation: {len(sf)} genes (symbols)")


# ── Query mygene.info ────────────────────────────────────────────────────────

mg = mygene.MyGeneInfo()

def query_by_ensembl(ensembl_ids: list, batch_size: int = 1000) -> dict:
    """Query mygene by Ensembl gene IDs."""
    results = {}
    for i in range(0, len(ensembl_ids), batch_size):
        batch = ensembl_ids[i:i+batch_size]
        out = mg.querymany(batch, scopes="ensembl.gene",
                           fields="symbol,type_of_gene,uniprot.Swiss-Prot,"
                                  "refseq.protein,ensembl.protein",
                           species="human", returnall=True)
        for r in out["out"]:
            results[r["query"]] = r
        if i + batch_size < len(ensembl_ids):
            time.sleep(0.5)
    return results


def query_by_symbol(symbols: list, batch_size: int = 1000) -> dict:
    """Query mygene by gene symbols."""
    results = {}
    for i in range(0, len(symbols), batch_size):
        batch = symbols[i:i+batch_size]
        out = mg.querymany(batch, scopes="symbol",
                           fields="ensembl.gene,type_of_gene,uniprot.Swiss-Prot,"
                                  "refseq.protein,ensembl.protein",
                           species="human", returnall=True)
        for r in out["out"]:
            if r["query"] not in results:  # take first hit per symbol
                results[r["query"]] = r
        if i + batch_size < len(symbols):
            time.sleep(0.5)
    return results


# Check for cached results
cache_file = CACHE / "E1_mygene_cache.json"
if cache_file.exists():
    print("\nLoading cached mygene results...")
    with open(cache_file) as f:
        cached = json.load(f)
    gf_results = cached["geneformer"]
    sg_results = cached["scgpt"]
    sf_results = cached["scfoundation"]
else:
    print("\nQuerying mygene.info (this may take a minute)...")

    print("  Geneformer (by Ensembl ID)...")
    gf_results = query_by_ensembl(gf["ensembl_id"].tolist())

    print("  scGPT (by symbol)...")
    sg_results = query_by_symbol(sg["gene"].tolist())

    print("  scFoundation (by symbol)...")
    sf_results = query_by_symbol(sf["gene"].tolist())

    # Cache results
    with open(cache_file, "w") as f:
        json.dump({
            "geneformer": gf_results,
            "scgpt": sg_results,
            "scfoundation": sf_results,
        }, f)
    print(f"  Cached to {cache_file}")


# ── Build feasibility table ──────────────────────────────────────────────────

def extract_mapping(result: dict) -> dict:
    """Extract protein mapping info from a mygene result."""
    if "notfound" in result and result["notfound"]:
        return {"found": False, "biotype": None, "has_protein": False,
                "uniprot_id": None, "n_proteins": 0}

    biotype = result.get("type_of_gene", None)

    # Check for protein IDs
    uniprot = result.get("uniprot", {})
    swissprot = uniprot.get("Swiss-Prot", None) if isinstance(uniprot, dict) else None
    if isinstance(swissprot, list):
        swissprot = swissprot[0]  # take canonical

    ensembl_proteins = result.get("ensembl", {})
    if isinstance(ensembl_proteins, dict):
        ensembl_proteins = ensembl_proteins.get("protein", [])
    elif isinstance(ensembl_proteins, list):
        # Multiple Ensembl entries — flatten proteins
        prots = []
        for e in ensembl_proteins:
            if isinstance(e, dict):
                p = e.get("protein", [])
                if isinstance(p, str):
                    prots.append(p)
                elif isinstance(p, list):
                    prots.extend(p)
        ensembl_proteins = prots
    if isinstance(ensembl_proteins, str):
        ensembl_proteins = [ensembl_proteins]

    refseq_proteins = result.get("refseq", {})
    if isinstance(refseq_proteins, dict):
        refseq_proteins = refseq_proteins.get("protein", [])
    if isinstance(refseq_proteins, str):
        refseq_proteins = [refseq_proteins]

    n_proteins = len(set(
        ([swissprot] if swissprot else []) +
        (ensembl_proteins or []) +
        (refseq_proteins or [])
    ))

    has_protein = swissprot is not None or bool(ensembl_proteins) or bool(refseq_proteins)

    return {
        "found": True,
        "biotype": biotype,
        "has_protein": has_protein,
        "uniprot_id": swissprot,
        "n_proteins": n_proteins,
    }


print("\nBuilding feasibility table...")

rows = []

# Geneformer
for _, row in gf.iterrows():
    gene = row["gene"]
    ens_id = row["ensembl_id"]
    is_outlier = row["is_outlier"]
    result = gf_results.get(ens_id, {"notfound": True})
    mapping = extract_mapping(result)

    rows.append({
        "model": "Geneformer",
        "vocab_id": ens_id,
        "gene_symbol": gene,
        "ensembl_gene": ens_id,
        "biotype": mapping["biotype"],
        "has_protein": mapping["has_protein"],
        "uniprot_id": mapping["uniprot_id"],
        "in_outlier_set": bool(is_outlier),
        "in_robust_core": gene in robust_cores.get("Geneformer", []),
        "gene_class": assign_gene_class(gene),
        "mapping_status": "mapped" if mapping["has_protein"] else
                         ("found_no_protein" if mapping["found"] else "not_found"),
    })

# scGPT
for _, row in sg.iterrows():
    gene = row["gene"]
    is_outlier = row["is_outlier"]
    result = sg_results.get(gene, {"notfound": True})
    mapping = extract_mapping(result)

    ens_gene = None
    if mapping["found"]:
        ens = result.get("ensembl", {})
        if isinstance(ens, dict):
            ens_gene = ens.get("gene", None)
        elif isinstance(ens, list) and ens:
            ens_gene = ens[0].get("gene", None) if isinstance(ens[0], dict) else None

    rows.append({
        "model": "scGPT",
        "vocab_id": gene,
        "gene_symbol": gene,
        "ensembl_gene": ens_gene,
        "biotype": mapping["biotype"],
        "has_protein": mapping["has_protein"],
        "uniprot_id": mapping["uniprot_id"],
        "in_outlier_set": bool(is_outlier),
        "in_robust_core": gene in robust_cores.get("scGPT", []),
        "gene_class": assign_gene_class(gene),
        "mapping_status": "mapped" if mapping["has_protein"] else
                         ("found_no_protein" if mapping["found"] else "not_found"),
    })

# scFoundation
for _, row in sf.iterrows():
    gene = row["gene"]
    is_outlier = row["is_outlier"]
    result = sf_results.get(gene, {"notfound": True})
    mapping = extract_mapping(result)

    ens_gene = None
    if mapping["found"]:
        ens = result.get("ensembl", {})
        if isinstance(ens, dict):
            ens_gene = ens.get("gene", None)
        elif isinstance(ens, list) and ens:
            ens_gene = ens[0].get("gene", None) if isinstance(ens[0], dict) else None

    rows.append({
        "model": "scFoundation",
        "vocab_id": gene,
        "gene_symbol": gene,
        "ensembl_gene": ens_gene,
        "biotype": mapping["biotype"],
        "has_protein": mapping["has_protein"],
        "uniprot_id": mapping["uniprot_id"],
        "in_outlier_set": bool(is_outlier),
        "in_robust_core": gene in robust_cores.get("scFoundation", []),
        "gene_class": assign_gene_class(gene),
        "mapping_status": "mapped" if mapping["has_protein"] else
                         ("found_no_protein" if mapping["found"] else "not_found"),
    })

feasibility = pd.DataFrame(rows)
feasibility.to_csv(OUT / "E1_feasibility_table.csv", index=False)
print(f"  Feasibility table: {len(feasibility)} rows")


# ── Bias audit ───────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("  E1 STAGE 1 — MAPPING FEASIBILITY & BIAS AUDIT")
print("=" * 70)

audit_rows = []

for model_name in ["Geneformer", "scGPT", "scFoundation"]:
    mdf = feasibility[feasibility["model"] == model_name].copy()

    # Classify as protein-coding or not
    mdf["is_coding"] = mdf["biotype"] == "protein-coding"

    n_total = len(mdf)
    n_coding = mdf["is_coding"].sum()
    n_outlier = mdf["in_outlier_set"].sum()
    n_outlier_coding = (mdf["in_outlier_set"] & mdf["is_coding"]).sum()
    n_outlier_noncoding = (mdf["in_outlier_set"] & ~mdf["is_coding"]).sum()

    # Mapping rates
    n_mapped = (mdf["mapping_status"] == "mapped").sum()
    n_mapped_coding = (mdf["is_coding"] & (mdf["mapping_status"] == "mapped")).sum()
    n_outlier_mapped = (mdf["in_outlier_set"] & (mdf["mapping_status"] == "mapped")).sum()
    n_outlier_coding_mapped = (mdf["in_outlier_set"] & mdf["is_coding"] &
                                (mdf["mapping_status"] == "mapped")).sum()

    # Mapping rate for protein-coding outliers (the gate metric)
    if n_outlier_coding > 0:
        outlier_coding_map_rate = n_outlier_coding_mapped / n_outlier_coding
    else:
        outlier_coding_map_rate = 0

    # Fisher's exact: coding/non-coding × outlier/non-outlier
    a = int((mdf["in_outlier_set"] & ~mdf["is_coding"]).sum())   # outlier & non-coding
    b = int((mdf["in_outlier_set"] & mdf["is_coding"]).sum())    # outlier & coding
    c = int((~mdf["in_outlier_set"] & ~mdf["is_coding"]).sum())  # non-outlier & non-coding
    d = int((~mdf["in_outlier_set"] & mdf["is_coding"]).sum())   # non-outlier & coding

    # OR for non-coding enrichment among outliers
    if b > 0 and c > 0:
        fisher_or, fisher_p = stats.fisher_exact([[a, b], [c, d]],
                                                   alternative="two-sided")
    else:
        fisher_or, fisher_p = np.nan, np.nan

    print(f"\n{'='*70}")
    print(f"  {model_name}")
    print(f"{'='*70}")
    print(f"  Total genes:          {n_total}")
    print(f"  Protein-coding:       {n_coding} ({100*n_coding/n_total:.1f}%)")
    print(f"  Outliers:             {n_outlier}")
    print(f"    coding:             {n_outlier_coding}")
    print(f"    non-coding/unknown: {n_outlier_noncoding}")
    print(f"  Overall mapping rate: {n_mapped}/{n_total} ({100*n_mapped/n_total:.1f}%)")
    print(f"  Coding mapping rate:  {n_mapped_coding}/{n_coding} "
          f"({100*n_mapped_coding/n_coding:.1f}%)" if n_coding > 0 else "  N/A")
    print(f"  Outlier mapping rate: {n_outlier_mapped}/{n_outlier} "
          f"({100*n_outlier_mapped/n_outlier:.1f}%)")
    print(f"  Outlier CODING map:   {n_outlier_coding_mapped}/{n_outlier_coding} "
          f"({100*outlier_coding_map_rate:.1f}%)" if n_outlier_coding > 0 else "  N/A")
    print(f"\n  Bias check (non-coding enrichment among outliers):")
    print(f"    2×2: outlier×non-coding={a}, outlier×coding={b}, "
          f"non-outlier×non-coding={c}, non-outlier×coding={d}")
    print(f"    Fisher OR={fisher_or:.2f}, p={fisher_p:.2e}")

    # Per gene-class breakdown
    print(f"\n  Per gene-class mapping rates:")
    for gc in ["ribosomal", "mitochondrial", "constrained", "disease", "other"]:
        gc_mask = mdf["gene_class"] == gc
        gc_n = gc_mask.sum()
        gc_outlier = (gc_mask & mdf["in_outlier_set"]).sum()
        gc_mapped = (gc_mask & (mdf["mapping_status"] == "mapped")).sum()
        if gc_n > 0:
            print(f"    {gc:15s}  n={gc_n:6d}  outlier={gc_outlier:4d}  "
                  f"mapped={gc_mapped:5d} ({100*gc_mapped/gc_n:.0f}%)")

    # Gate assessment
    gate_map_rate = outlier_coding_map_rate
    gate_or = fisher_or

    if gate_map_rate >= 0.75 and gate_or < 3.0:
        gate = "PROCEED"
    elif gate_map_rate >= 0.50 and gate_or < 5.0:
        gate = "SCOPE-RESTRICT"
    else:
        gate = "STOP"

    print(f"\n  GATE: {gate}")
    print(f"    Protein-coding outlier mapping: {100*gate_map_rate:.1f}% "
          f"(threshold: ≥75%)")
    print(f"    Non-coding enrichment OR: {gate_or:.2f} "
          f"(threshold: <3.0)")

    audit_rows.append({
        "model": model_name,
        "n_total": n_total,
        "n_coding": n_coding,
        "n_outlier": n_outlier,
        "n_outlier_coding": n_outlier_coding,
        "n_outlier_noncoding": n_outlier_noncoding,
        "overall_mapping_rate": n_mapped / n_total,
        "coding_mapping_rate": n_mapped_coding / n_coding if n_coding > 0 else np.nan,
        "outlier_mapping_rate": n_outlier_mapped / n_outlier if n_outlier > 0 else np.nan,
        "outlier_coding_mapping_rate": outlier_coding_map_rate,
        "fisher_OR_noncoding": fisher_or,
        "fisher_p_noncoding": fisher_p,
        "gate": gate,
    })


# ── Overall gate ─────────────────────────────────────────────────────────────

audit_df = pd.DataFrame(audit_rows)
audit_df.to_csv(OUT / "E1_bias_audit.csv", index=False)

gates = audit_df["gate"].tolist()
if all(g == "PROCEED" for g in gates):
    overall = "PROCEED"
elif any(g == "STOP" for g in gates):
    overall = "STOP (at least one model fails)"
else:
    overall = "SCOPE-RESTRICT"

print(f"\n\n{'='*70}")
print(f"  OVERALL E1 STAGE 1 GATE: {overall}")
print(f"{'='*70}")
for _, r in audit_df.iterrows():
    print(f"  {r['model']:15s}  {r['gate']:16s}  "
          f"coding-outlier-map={100*r['outlier_coding_mapping_rate']:.0f}%  "
          f"OR={r['fisher_OR_noncoding']:.2f}")

# Save gate verdict
gate_verdict = {
    "overall": overall,
    "panel_provenance": panel_provenance(),
    "per_model": {r["model"]: {
        "gate": r["gate"],
        "outlier_coding_mapping_rate": r["outlier_coding_mapping_rate"],
        "fisher_OR_noncoding": r["fisher_OR_noncoding"],
    } for _, r in audit_df.iterrows()},
}
with open(OUT / "E1_stage1_gate.json", "w") as f:
    json.dump(gate_verdict, f, indent=2)

print(f"\n  Saved to {OUT}/:")
print(f"    E1_feasibility_table.csv ({len(feasibility)} rows)")
print(f"    E1_bias_audit.csv")
print(f"    E1_stage1_gate.json")
print(f"  Cached to {CACHE}/:")
print(f"    E1_mygene_cache.json")
