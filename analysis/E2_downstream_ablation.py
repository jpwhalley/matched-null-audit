"""
E2 — Expression-matched downstream ablation.

Tests whether removing geometric-outlier genes from Geneformer's input
degrades cell-type annotation more than removing expression-matched,
class-stratified control genes.

Analysis parameters, fixed before any ablation results were inspected:
  - Model: Geneformer V2-104M (most stable outlier set per E3)
  - top-k: 50 (k=25 and k=100 were also run; neither is reported in the
    manuscript, see MANUSCRIPT_TRACEABILITY.md)
  - Data: PBMC3k + Tabula Sapiens immune subset
  - Treatment: delete gene tokens from cell sequences
  - Controls: class-stratified matching on log1p(expression), breadth,
    log1p(gene_length) with pool-level z-scored distances.
    200 matched-control draws per treatment set.
  - Primary metric: macro-averaged F1 (retrained linear probe, 5-fold CV)
  - Fixed-probe metric: train on baseline, predict on ablated (no recalibration)
  - Clustering metrics: k-means ARI and NMI vs true labels
  - Sensitivity: drop ribosomal + mitochondrial outliers, re-test against
    *own matched null* (36-gene treatment → 36-gene controls)

Usage:
  python E2_downstream_ablation.py --setup        # download data, tokenize,
                                                  # build matched controls
  python E2_downstream_ablation.py --baseline     # baseline embeddings + metrics
  python E2_downstream_ablation.py --ablation     # treatment + 2×200 controls
  python E2_downstream_ablation.py --evaluate     # compare, write verdict
  python E2_downstream_ablation.py --all          # everything

Requirements (already in pyproject.toml except cellxgene-census):
  torch, transformers, anndata, scanpy, scikit-learn, scipy, pandas

Outputs (all in outputs/):
  E2_baseline_{dataset}.json           — baseline metrics
  E2_treatment_genes.csv               — treatment gene set with classes
  E2_ablation_{dataset}.json           — full ablation results
  E2_verdict_{dataset}.json            — structured gate verdict
  E2_matched_controls_*_balance.csv    — matching diagnostics (SMD table)
"""

import argparse
import collections
import functools
import hashlib
import pathlib
import json
import os
import pickle
import re
import sys
import warnings
from pathlib import Path

# ── Environment pinning (MUST run before numpy/torch import) ────────────────
# H1: the PBMC3k baseline macro-F1 differed by 0.00058 between machines while
# the embeddings were bit-identical, i.e. the probe fit is environment
# sensitive. BLAS thread count is the most likely lever, so pin it by default
# and record the full fingerprint alongside every result. Set E2_PIN_THREADS=0
# to opt out (not recommended for numbers destined for a manuscript).
_PIN_THREADS = os.environ.get("E2_PIN_THREADS", "1")
if _PIN_THREADS != "0":
    for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[_v] = _PIN_THREADS

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.cluster import KMeans
from sklearn.metrics import (f1_score, adjusted_rand_score,
                             normalized_mutual_info_score)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder

# Repository-relative paths: this script reads and writes only inside the
# repository, under data/, outputs/ and cache/. The Geneformer checkpoint is
# expected at BASE/"Geneformer" and is not shipped, which is why a fresh
# clone cannot run the ablation without acquiring it first (DATA_MANIFEST.md).
# Repository-relative paths. Scripts live in analysis/; everything they read
# and write is inside this repository.
REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUT = REPO / "outputs"
CACHE = REPO / "cache"
for _d in (OUT, CACHE):
    _d.mkdir(parents=True, exist_ok=True)
BASE = REPO  # legacy alias


try:
    import torch
except ImportError:
    torch = None  # torch only needed for Steps 2-3

warnings.filterwarnings("ignore")

OUT  = OUT

# ── Environment fingerprint ─────────────────────────────────────────────────

def environment_fingerprint(device=None):
    """Record everything that could move a probe fit by ~1e-4 (see H1).

    Saved alongside every baseline and ablation result so a number can be
    traced to the environment that produced it. If two runs disagree, diff
    these blocks first.
    """
    import platform
    fp = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "threads_pinned": _PIN_THREADS,
        "thread_env": {v: os.environ.get(v) for v in
                       ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                        "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                        "NUMEXPR_NUM_THREADS")},
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    try:
        import scipy
        fp["scipy"] = scipy.__version__
    except Exception:
        pass
    try:
        import sklearn
        fp["sklearn"] = sklearn.__version__
    except Exception:
        pass
    try:
        import torch as _t
        fp["torch"] = _t.__version__
        fp["torch_num_threads"] = _t.get_num_threads()
    except Exception:
        pass
    try:
        import threadpoolctl
        fp["threadpools"] = [
            {k: str(d.get(k)) for k in
             ("user_api", "internal_api", "prefix", "version",
              "num_threads")}
            for d in threadpoolctl.threadpool_info()]
    except Exception:
        fp["threadpools"] = "threadpoolctl not installed (pip install "\
                            "threadpoolctl for full BLAS detail)"
    if device is not None:
        fp["device"] = str(device)
    return fp


# ── Locked parameters ───────────────────────────────────────────────────────

TOP_K = 50
N_BOOTSTRAP = 200
SUBSAMPLE_N = None  # set by --subsample; stratified by cell type
MODEL_DIR = BASE / "Geneformer" / "Geneformer-V2-104M"
FORWARD_BATCH_SIZE = 8  # max cells per batch (see ATTENTION_SQ_BUDGET)
# Peak CPU attention memory scales with batch x seq_len^2, and dynamic padding
# means a batch costs as much as its LONGEST member. Cell-count batching alone
# therefore spikes wherever long cells cluster (the Smart-seq2 block in Tabula
# Sapiens), which killed the 20k-cell run at a reproducible index (~15.2k).
# Cap batch x max_len^2 instead: batch falls to 1 cell at 4096 tokens.
# Lower this if an OOM recurs; raise it for speed if memory headroom allows.
ATTENTION_SQ_BUDGET = 16 * 1024 ** 2
N_CV_FOLDS = 5
# Cell types below this are dropped from Tabula Sapiens BEFORE analysis.
# macro-F1 weights every class equally, so a long tail of 5-9 cell classes
# would dominate the metric's variance while holding ~2% of the data.
MIN_CELLS_PER_TYPE = 10
RANDOM_SEED = 42

# ── Gene class assignment (same as E1/E3) ───────────────────────────────────

# Full-coverage annotation table, checksum-pinned in data/CHECKSUMS.json. The
# length-matched design needs near-complete gene lengths: the coverage assert
# below refuses any table below MIN_LENGTH_COVERAGE, because median-imputing a
# large fraction of lengths silently corrupts the balance statistics.
TABLE_S1_PATH = DATA / "Table_S1.csv"
MIN_LENGTH_COVERAGE = 0.99

table_s1 = pd.read_csv(TABLE_S1_PATH)

_len_cov = float(table_s1["gene_length_bp"].notna().mean())
if _len_cov < MIN_LENGTH_COVERAGE:
    raise RuntimeError(
        f"Gene-length coverage in {TABLE_S1_PATH} is {_len_cov:.4f}, below the "
        f"required {MIN_LENGTH_COVERAGE}. Matching on a median-imputed length "
        f"silently breaks the balance the design depends on. Point "
        f"TABLE_S1_PATH at a full-coverage table.")


def _table_s1_provenance():
    """Hash + coverage, recorded in every output so a run is traceable."""
    import hashlib
    h = hashlib.sha256()
    with open(TABLE_S1_PATH, "rb") as _f:
        for _chunk in iter(lambda: _f.read(1 << 20), b""):
            h.update(_chunk)
    return {"path": str(TABLE_S1_PATH.relative_to(REPO)),
            "sha256_16": h.hexdigest()[:16],
            "n_genes": int(len(table_s1)),
            "gene_length_coverage": round(_len_cov, 4)}

# Pre-compute constrained and disease gene sets (avoid recomputing per call)
_CONSTRAINED_GENES = set(
    table_s1.loc[(table_s1["pLI"] > 0.9) | (table_s1["LOEUF"] < 0.35),
                 "gene_symbol"].str.upper()
)
_DISEASE_GENES = set(
    table_s1.loc[table_s1["clinvar_disease"] == True,
                 "gene_symbol"].str.upper()
)


# Ribosomal panel: pinned HGNC gene groups 728, 729 and 646, resolved by
# Ensembl gene ID rather than by symbol, since symbol matching admits RPS6K
# signalling kinases and misses members carrying obsolete symbols. See
# analysis/_ribosomal_panel.py and data/ribosomal_panel_provenance.json.
# Every matched-control cache records the panel SHA-256 and is refused if it
# differs.
from _ribosomal_panel import (ribosomal_symbols, panel_provenance,  # noqa: E402
                              PANEL_SHA256 as RIBOSOMAL_PANEL_SHA256)
_RIBOSOMAL = ribosomal_symbols(table_s1)


def assign_gene_class(sym: str) -> str:
    """Class labels for the matched-control strata (mutually exclusive)."""
    sym_u = sym.upper()
    if sym_u.startswith("MT-"):
        return "mitochondrial"
    if sym_u in _RIBOSOMAL:
        return "ribosomal"
    if sym_u in _CONSTRAINED_GENES:
        return "constrained"
    if sym_u in _DISEASE_GENES:
        return "disease"
    return "other"


# ── Helpers ──────────────────────────────────────────────────────────────────

def load_geneformer_dicts():
    """Load token dictionary and gene median dictionary."""
    gf_pkg = BASE / "Geneformer" / "geneformer"
    with open(gf_pkg / "token_dictionary_gc104M.pkl", "rb") as f:
        token_dict = pickle.load(f)  # ensembl_id -> token_id
    with open(gf_pkg / "gene_median_dictionary_gc104M.pkl", "rb") as f:
        gene_median_dict = pickle.load(f)
    return token_dict, gene_median_dict


def get_treatment_genes(k=TOP_K):
    """Return top-k outlier genes with their Ensembl IDs and token IDs."""
    gf_geom = pd.read_csv(DATA / "gene_embedding_geometry.csv")
    gf_geom = gf_geom[~gf_geom["gene"].isin(["<pad>", "<mask>", "<cls>", "<eos>"])]
    sc_col = ("anomaly_score_with_isolation"
              if "anomaly_score_with_isolation" in gf_geom.columns
              else "anomaly_score")
    top_k = gf_geom.nlargest(k, sc_col).copy()
    top_k["gene_class"] = top_k["gene"].map(assign_gene_class)
    return top_k


def tokenize_cell_from_adata(gene_vector, gene_ensembl_ids,
                             token_dict, gene_median_dict,
                             special_token=True):
    """Tokenize a single cell: rank-value encoding with optional CLS/EOS.

    gene_vector: raw count vector (1D array)
    gene_ensembl_ids: matching Ensembl IDs
    token_dict: ensembl_id -> token_id
    gene_median_dict: ensembl_id -> median expression
    """
    # Filter to genes in vocabulary
    mask = np.array([eid in token_dict and eid in gene_median_dict
                     for eid in gene_ensembl_ids])
    counts = gene_vector[mask]
    eids = gene_ensembl_ids[mask]

    # Normalise to 10k and divide by median
    total = counts.sum()
    if total == 0:
        return np.array([], dtype=np.int64)
    norm = (counts / total) * 10000
    median_scaled = np.array([
        norm[i] / gene_median_dict[eids[i]]
        for i in range(len(eids))
        if norm[i] > 0 and gene_median_dict.get(eids[i], 0) > 0
    ])
    valid_eids = np.array([
        eids[i] for i in range(len(eids))
        if norm[i] > 0 and gene_median_dict.get(eids[i], 0) > 0
    ])

    if len(median_scaled) == 0:
        return np.array([], dtype=np.int64)

    # Rank by median-scaled value (descending)
    sorted_idx = np.argsort(-median_scaled)
    token_ids = np.array([token_dict[valid_eids[i]] for i in sorted_idx],
                         dtype=np.int64)

    # Truncate to model input size - 2 (for CLS + EOS)
    max_len = 4094 if special_token else 4096
    token_ids = token_ids[:max_len]

    if special_token:
        cls_id = token_dict.get("<cls>", 2)
        eos_id = token_dict.get("<eos>", 3)
        token_ids = np.concatenate([[cls_id], token_ids, [eos_id]])

    return token_ids


def delete_gene_tokens(token_seq, tokens_to_delete, pad_id=0):
    """Remove specified tokens from sequence, pad to maintain length."""
    delete_set = set(tokens_to_delete)
    kept = [t for t in token_seq if t not in delete_set]
    n_pad = len(token_seq) - len(kept)
    return np.array(kept + [pad_id] * n_pad, dtype=np.int64)


# ── Step 1: Setup ────────────────────────────────────────────────────────────

def _strip_ensembl_version(ids):
    """ENSG00000000003.14 -> ENSG00000000003"""
    return np.array([eid.split(".")[0] if isinstance(eid, str) else eid
                     for eid in ids])


def _ts_celltype_col(adata):
    """First available cell-type annotation column, or None."""
    for candidate in ["cell_type", "cell_ontology_class", "free_annotation",
                      "compartment", "louvain"]:
        if candidate in adata.obs.columns:
            return candidate
    return None


def _stratified_subsample(adata, n_target, ct_col, seed=RANDOM_SEED):
    """Subsample to ~n_target cells preserving cell-type proportions.

    Proportional allocation with a floor of 1 cell per type, so rare types
    are not silently dropped -- losing a class entirely would change the
    classification task rather than shrink it.
    """
    labels = pd.Series(adata.obs[ct_col].values).astype(str)
    n_total = len(labels)
    if n_target >= n_total:
        print(f"    Subsample target {n_target} >= {n_total} cells; "
              f"using all cells.")
        return adata

    rng = np.random.default_rng(seed)
    counts = labels.value_counts()

    # Purely proportional. No floor is applied: classes that fall below
    # MIN_CELLS_PER_TYPE are dropped by the downstream filter anyway, and
    # inflating them here would only create classes that are then discarded.
    floor = 1
    take = {}
    for ct, n in counts.items():
        take[ct] = max(1, int(round(n * n_target / n_total)))

    # Trim/pad the largest types so the total lands on n_target
    total = sum(take.values())
    order = list(counts.index)
    i = 0
    while total != n_target and i < 100000:
        ct = order[i % len(order)]
        if total > n_target and take[ct] > min(floor, int(counts[ct])):
            take[ct] -= 1; total -= 1
        elif total < n_target and take[ct] < counts[ct]:
            take[ct] += 1; total += 1
        i += 1

    idx = []
    for ct, k in take.items():
        pool = np.flatnonzero((labels == ct).values)
        idx.extend(rng.choice(pool, size=min(k, len(pool)), replace=False))
    idx = np.sort(np.array(idx))

    sub = adata[idx].copy()
    n_types_before, n_types_after = counts.size, \
        pd.Series(sub.obs[ct_col].values).astype(str).nunique()
    print(f"    Stratified subsample: {n_total} -> {sub.shape[0]} cells, "
          f"{n_types_before} -> {n_types_after} cell types (seed={seed})")
    if n_types_after < n_types_before:
        print(f"    WARNING: {n_types_before - n_types_after} cell type(s) "
              f"lost despite the floor of {floor}.")
    per_class = pd.Series(sub.obs[ct_col].values).astype(str).value_counts()
    surviving = per_class[per_class >= MIN_CELLS_PER_TYPE]
    print(f"    Smallest class after subsample: {per_class.min()} cells")
    print(f"    Downstream filter (>={MIN_CELLS_PER_TYPE} cells/type) will "
          f"retain {len(surviving)}/{len(per_class)} classes, "
          f"{surviving.sum()}/{sub.shape[0]} cells "
          f"({100 * surviving.sum() / sub.shape[0]:.0f}%)")
    print(f"    NOTE: macro-F1 weights classes equally, so the dropped "
          f"classes would have contributed "
          f"{100 * (len(per_class) - len(surviving)) / len(per_class):.0f}% "
          f"of the metric on "
          f"{100 * (sub.shape[0] - surviving.sum()) / sub.shape[0]:.0f}% of "
          f"the cells. Report Tabula Sapiens as a "
          f"{len(surviving)}-class filtered immune subset.")
    print(f"    The TS null band will still be wider than PBMC3k's. A wider "
          f"band favours the 'inside' verdict, so treat a TS null as "
          f"supportive only, never as independent confirmation.")
    return sub


def _tokenize_and_cache_dataset(adata, ensembl_ids, ct_col,
                                 token_dict, gene_median_dict,
                                 dataset_name):
    """Tokenize an AnnData, compute gene stats, save caches. Returns gene_stats_df."""
    import scipy.sparse as sp_sparse

    # Strip Ensembl version suffixes (CellxGene uses ENSG*.14, GF uses ENSG*)
    ensembl_ids = _strip_ensembl_version(ensembl_ids)

    # Use raw counts when available (CellxGene .X is often normalised)
    if hasattr(adata, "raw") and adata.raw is not None:
        print("    Using .raw.X (raw counts)")
        X = adata.raw.X
        # raw may have different var; re-derive ensembl IDs
        raw_var = adata.raw.var
        for col in ["feature_id", "ensembl_id"]:
            if col in raw_var.columns:
                ensembl_ids = _strip_ensembl_version(raw_var[col].values)
                break
        else:
            if raw_var.index[0].startswith("ENSG"):
                ensembl_ids = _strip_ensembl_version(raw_var.index.values)
    else:
        X = adata.X

    if sp_sparse.issparse(X):
        X = X.toarray()

    # Sanity: check vocabulary overlap before tokenizing 20k cells
    n_in_vocab = sum(1 for eid in ensembl_ids
                     if eid in token_dict and eid in gene_median_dict)
    print(f"    Genes in Geneformer vocabulary: {n_in_vocab}/{len(ensembl_ids)}")
    if n_in_vocab < 100:
        print("    WARNING: Very few genes match vocabulary.")
        print(f"    Sample IDs: {list(ensembl_ids[:5])}")
        print(f"    Sample dict keys: {list(token_dict.keys())[:5]}")

    print(f"    Tokenizing {X.shape[0]} cells...")
    tokenized = []
    for i in range(X.shape[0]):
        tokens = tokenize_cell_from_adata(
            X[i], ensembl_ids, token_dict, gene_median_dict
        )
        tokenized.append(tokens)
        if (i + 1) % 5000 == 0:
            print(f"      {i+1}/{X.shape[0]}")

    # Diagnostic: token length distribution
    lengths = [len(t) for t in tokenized]
    med_len = np.median(lengths)
    print(f"    Token lengths: median={med_len:.0f}, "
          f"min={min(lengths)}, max={max(lengths)}")
    if med_len < 10:
        raise RuntimeError(
            f"Tokenization failed for {dataset_name}: median token "
            f"length = {med_len:.0f}. Check that Ensembl IDs match "
            f"Geneformer vocabulary (no version suffixes) and that "
            f".X contains raw counts (not normalised values)."
        )

    print("    Computing gene expression statistics...")
    gene_stats = []
    for j, eid in enumerate(ensembl_ids):
        if eid not in token_dict:
            continue
        col = X[:, j]
        gene_stats.append({
            "ensembl_id": eid,
            "token_id": token_dict[eid],
            "expr_mean": float(col.mean()),
            "expr_breadth": float((col > 0).mean()),
        })
    gene_stats_df = pd.DataFrame(gene_stats)

    cache_data = {
        "tokenized": [t.tolist() for t in tokenized],
        "cell_types": adata.obs[ct_col].values.tolist(),
        "n_cells": len(tokenized),
        "subsample_n": SUBSAMPLE_N if SUBSAMPLE_N else "unset",
    }
    with open(CACHE / f"E2_{dataset_name}_tokenized.json", "w") as f:
        json.dump(cache_data, f)
    gene_stats_df.to_csv(CACHE / f"E2_{dataset_name}_gene_stats.csv",
                         index=False)

    print(f"    Saved: {len(tokenized)} cells tokenized")
    print(f"    Median token length: "
          f"{np.median([len(t) for t in tokenized]):.0f}")
    return gene_stats_df


def _download_tabula_sapiens():
    """Download Tabula Sapiens - Immune from CellxGene. Returns path or None."""
    import anndata as ad

    # CellxGene dataset: "Tabula Sapiens - Immune"
    TS_DATASET_ID = "78b60b70-129a-4a6d-b15f-825b241eec66"

    # Check for user-placed file first
    manual_paths = [
        CACHE / "tabula_sapiens_immune.h5ad",
        CACHE / "E2_tabula_sapiens_immune.h5ad",
        DATA / "tabula_sapiens_immune.h5ad",
        DATA / "tabula_sapiens.h5ad",
    ]
    for mp in manual_paths:
        if mp.exists():
            print(f"    Found local file: {mp.name}")
            return mp

    # Method 1: cellxgene-census (streams only rows needed)
    try:
        import cellxgene_census
        print("    Downloading via cellxgene-census (streaming)...")
        with cellxgene_census.open_soma() as census:
            ts = cellxgene_census.get_anndata(
                census,
                organism="Homo sapiens",
                obs_value_filter=(
                    f"dataset_id == '{TS_DATASET_ID}'"
                ),
            )
        out_path = CACHE / "tabula_sapiens_immune_full.h5ad"
        ts.write(out_path)
        print(f"    Downloaded {ts.shape[0]} cells via census API")
        return out_path
    except ImportError:
        print("    cellxgene-census not installed (try: uv add cellxgene-census)")
    except Exception as e:
        print(f"    cellxgene-census failed: {e}")

    # Method 2: direct download via CellxGene API (19.8 GB)
    print("    Trying direct download from CellxGene (~19.8 GB)...")
    print("    This will take a while on a typical connection.")
    import urllib.request

    # Get the download URL via the assets API
    try:
        assets_url = (f"https://api.cellxgene.cziscience.com/dp/v1/"
                      f"datasets/{TS_DATASET_ID}/assets")
        assets = json.loads(urllib.request.urlopen(assets_url, timeout=30).read())
        h5ad_asset = next(
            (a for a in assets.get("assets", [])
             if a.get("filetype") == "H5AD"), None
        )
        if h5ad_asset is None:
            print("    ERROR: No H5AD asset found in CellxGene API response")
            return None

        # Get presigned download URL
        asset_info_url = (
            f"https://api.cellxgene.cziscience.com/dp/v1/"
            f"datasets/{TS_DATASET_ID}/asset/{h5ad_asset['id']}"
        )
        asset_info = json.loads(
            urllib.request.urlopen(asset_info_url, timeout=30).read()
        )
        dl_url = asset_info["url"]
        file_size = asset_info.get("file_size", 0)
        print(f"    Download URL obtained ({file_size / 1e9:.1f} GB)")

        dl_path = CACHE / "tabula_sapiens_immune_dl.h5ad"
        if dl_path.exists() and dl_path.stat().st_size > 1e9:
            print(f"    Partial download exists ({dl_path.stat().st_size / 1e9:.1f} GB)")
            # Check if it's complete
            if file_size > 0 and dl_path.stat().st_size >= file_size * 0.99:
                print("    Looks complete, using it.")
                return dl_path

        # Download with progress
        def _report(block_num, block_size, total_size):
            downloaded = block_num * block_size
            if total_size > 0:
                pct = downloaded / total_size * 100
                gb = downloaded / 1e9
                if block_num % 5000 == 0:
                    print(f"      {gb:.1f} GB / {total_size/1e9:.1f} GB "
                          f"({pct:.0f}%)")

        urllib.request.urlretrieve(dl_url, dl_path, reporthook=_report)
        print(f"    Download complete: {dl_path.stat().st_size / 1e9:.1f} GB")
        return dl_path

    except Exception as e:
        print(f"    Direct download failed: {e}")
        return None


def setup(datasets=None):
    import anndata as ad
    import scanpy as sc

    print("=" * 70)
    print("  STEP 1: Dataset setup")
    print("=" * 70)

    token_dict, gene_median_dict = load_geneformer_dicts()

    # ── Treatment gene set (needed for both datasets) ────────────────────
    print("\n  Defining treatment gene set (top-50)...")
    treatment = get_treatment_genes(TOP_K)
    treatment_info = treatment[["gene", "ensembl_id", "token_id",
                                "anomaly_score_with_isolation",
                                "gene_class"]].copy()
    treatment_info.to_csv(OUT / "E2_treatment_genes.csv", index=False)
    print(f"    Treatment genes: {len(treatment_info)}")
    print(f"    Class breakdown: "
          f"{dict(treatment_info['gene_class'].value_counts())}")

    # ── PBMC3k ───────────────────────────────────────────────────────────
    print("\n  ── PBMC3k ──")
    pbmc_tok_cache = CACHE / "E2_pbmc3k_tokenized.json"
    pbmc_stats_cache = CACHE / "E2_pbmc3k_gene_stats.csv"

    if pbmc_tok_cache.exists() and pbmc_stats_cache.exists():
        print("    Already tokenized (cached). Skipping.")
        pbmc_gene_stats = pd.read_csv(pbmc_stats_cache)
    else:
        pbmc_path = DATA / "pbmc3k_h5ad" / "pbmc3k.h5ad"
        adata = ad.read_h5ad(pbmc_path)
        print(f"    Shape: {adata.shape}")
        print(f"    Cell types: {dict(adata.obs['cell_type'].value_counts())}")

        if "ensembl_id" in adata.var.columns:
            ensembl_ids = adata.var["ensembl_id"].values
        elif adata.var.index[0].startswith("ENSG"):
            ensembl_ids = adata.var.index.values
        else:
            raise ValueError("Cannot find Ensembl IDs in PBMC3k")

        pbmc_gene_stats = _tokenize_and_cache_dataset(
            adata, ensembl_ids, "cell_type",
            token_dict, gene_median_dict, "pbmc3k"
        )

    # No-ribo/mito treatment subset (for sensitivity test)
    treatment_no_rm = treatment_info[
        ~treatment_info["gene_class"].isin(["ribosomal", "mitochondrial"])
    ].copy()
    print(f"    Sensitivity subset (no ribo/mito): {len(treatment_no_rm)} genes")

    # Build matched controls for PBMC3k — full and sensitivity
    pbmc_ctrl_cache = CACHE / "E2_matched_controls_pbmc3k.json"
    pbmc_sens_cache = CACHE / "E2_matched_controls_pbmc3k_no_ribo_mito.json"
    if pbmc_ctrl_cache.exists():
        _assert_control_spec(pbmc_ctrl_cache)
        print("    Matched controls (full) already built (cached).")
    else:
        print("    Building matched controls — full treatment (PBMC3k)...")
        build_matched_controls(treatment_info, pbmc_gene_stats, "pbmc3k")
    if pbmc_sens_cache.exists():
        _assert_control_spec(pbmc_sens_cache)
        print("    Matched controls (sensitivity) already built (cached).")
    else:
        print("    Building matched controls — no ribo/mito (PBMC3k)...")
        build_matched_controls(treatment_no_rm, pbmc_gene_stats, "pbmc3k",
                               label="no_ribo_mito")

    # ── Tabula Sapiens immune ────────────────────────────────────────────
    print("\n  ── Tabula Sapiens ──")
    if datasets is not None and 'tabula_sapiens' not in datasets:
        print('\n  Skipping Tabula Sapiens (not in --datasets)')
        print('\n  Setup complete.')
        return

    ts_tok_cache = CACHE / "E2_tabula_sapiens_tokenized.json"
    ts_stats_cache = CACHE / "E2_tabula_sapiens_gene_stats.csv"

    if ts_tok_cache.exists() and ts_stats_cache.exists():
        # Cache filenames are stable, so a cache built at a DIFFERENT
        # subsample size would be silently reused -- wrong gene statistics,
        # wrong matched controls, wrong cell count. Fail loudly instead.
        with open(ts_tok_cache) as _f:
            _cached_n = json.load(_f).get("subsample_n", "unset")
        _want = SUBSAMPLE_N if SUBSAMPLE_N else "unset"
        if _cached_n != _want:
            raise RuntimeError(
                f"Tabula Sapiens cache was built with subsample_n="
                f"{_cached_n!r} but this run requests {_want!r}.\n"
                f"Gene statistics and matched controls derived from it would "
                f"not correspond to the cells being analysed.\n"
                f"Delete the stale caches and re-run --setup:\n"
                f"  rm -f cache/E2_tabula_sapiens_tokenized.json \\\n"
                f"        cache/E2_tabula_sapiens_gene_stats.csv \\\n"
                f"        cache/E2_matched_controls_tabula_sapiens*.json\n"
                f"(Keep E2_tabula_sapiens_immune.h5ad - no need to re-download.)"
            )
        print(f"    Already tokenized (cached, subsample_n={_cached_n}). "
              f"Skipping.")
        ts_gene_stats = pd.read_csv(ts_stats_cache)
    else:
        ts_h5ad_cache = CACHE / "E2_tabula_sapiens_immune.h5ad"

        if ts_h5ad_cache.exists():
            print("    Loading cached h5ad...")
            ts = ad.read_h5ad(ts_h5ad_cache)
        else:
            ts_path = _download_tabula_sapiens()
            if ts_path is None:
                print("\n    *** Tabula Sapiens not available. ***")
                print("    Options:")
                print("      1. uv add cellxgene-census  (then re-run --setup)")
                print("      2. Download manually from CellxGene Discover:")
                print("         https://cellxgene.cziscience.com/collections/"
                      "e5f58829-1a66-40b5-a624-9046778e74f5")
                print("         Save as: cache/"
                      "tabula_sapiens_immune.h5ad")
                print("    PBMC3k setup is complete. Re-run --setup after "
                      "downloading.\n")
                return

            ts = ad.read_h5ad(ts_path)
            print(f"    Raw shape: {ts.shape}")

            # Subsample if very large (>20k cells)
            if ts.shape[0] > 20000:
                print(f"    Subsampling {ts.shape[0]} -> 20000 cells "
                      f"(seed={RANDOM_SEED})")
                np.random.seed(RANDOM_SEED)
                idx = np.random.choice(ts.shape[0], 20000, replace=False)
                ts = ts[idx].copy()

            ts.write(ts_h5ad_cache)
            print(f"    Saved subsampled h5ad: {ts.shape}")

        # --subsample: stratified reduction for tractability. Applied AFTER
        # the 20k cache so the expensive download/subsample is reused.
        if SUBSAMPLE_N:
            _ct = _ts_celltype_col(ts)
            if _ct is None:
                print("    WARNING: no cell-type column; cannot stratify. "
                      "Skipping subsample.")
            else:
                ts = _stratified_subsample(ts, SUBSAMPLE_N, _ct)

        # Identify cell type column
        ct_col = _ts_celltype_col(ts)
        if ct_col is None:
            print("    WARNING: No cell type column found. Skipping TS.")
            print("    Available columns:", list(ts.obs.columns[:20]))
            return

        print(f"    Cell type column: {ct_col}")
        cts = ts.obs[ct_col].value_counts()
        print(f"    Cell types: {len(cts)} unique")
        print(f"    Top 5: {dict(cts.head(5))}")

        # Filter to cell types with >= MIN_CELLS_PER_TYPE cells.
        # Applied AFTER any --subsample, so the threshold refers to the cells
        # actually analysed rather than to the full atlas.
        valid_cts = cts[cts >= MIN_CELLS_PER_TYPE].index
        ts = ts[ts.obs[ct_col].isin(valid_cts)].copy()
        _final = ts.obs[ct_col].value_counts()
        print(f"    After filtering (>={MIN_CELLS_PER_TYPE} cells/type): "
              f"{ts.shape[0]} cells, {len(valid_cts)} types "
              f"(smallest class {_final.min()})")
        if _final.min() < N_CV_FOLDS:
            raise RuntimeError(
                f"Smallest cell type has {_final.min()} cells but "
                f"{N_CV_FOLDS}-fold stratified CV needs >= {N_CV_FOLDS}. "
                f"Raise MIN_CELLS_PER_TYPE or the --subsample target.")

        # Get Ensembl IDs
        ts_ensembl = None
        for col in ["feature_id", "ensembl_id"]:
            if col in ts.var.columns:
                ts_ensembl = ts.var[col].values
                break
        if ts_ensembl is None and ts.var.index[0].startswith("ENSG"):
            ts_ensembl = ts.var.index.values
        if ts_ensembl is None:
            print("    WARNING: Cannot find Ensembl IDs in Tabula Sapiens")
            print("    Available var columns:", list(ts.var.columns[:20]))
            return

        ts_gene_stats = _tokenize_and_cache_dataset(
            ts, ts_ensembl, ct_col,
            token_dict, gene_median_dict, "tabula_sapiens"
        )

    # Build matched controls for TS — full and sensitivity
    ts_ctrl_cache = CACHE / "E2_matched_controls_tabula_sapiens.json"
    ts_sens_cache = CACHE / "E2_matched_controls_tabula_sapiens_no_ribo_mito.json"
    if ts_ctrl_cache.exists():
        _assert_control_spec(ts_ctrl_cache)
        with open(ts_ctrl_cache) as _f:
            _n_draws = len(json.load(_f))
        if _n_draws != N_BOOTSTRAP:
            raise RuntimeError(
                f"Cached Tabula Sapiens controls hold {_n_draws} draws but "
                f"N_BOOTSTRAP={N_BOOTSTRAP}. Delete "
                f"cache/E2_matched_controls_tabula_sapiens*.json "
                f"and re-run --setup.")
        print("    Matched controls (full) already built (cached).")
    else:
        print("    Building matched controls — full treatment (TS)...")
        build_matched_controls(treatment_info, ts_gene_stats,
                               "tabula_sapiens")
    if ts_sens_cache.exists():
        _assert_control_spec(ts_sens_cache)
        print("    Matched controls (sensitivity) already built (cached).")
    else:
        print("    Building matched controls — no ribo/mito (TS)...")
        build_matched_controls(treatment_no_rm, ts_gene_stats,
                               "tabula_sapiens", label="no_ribo_mito")

    print("\n  Setup complete.")


def _prepare_candidate_pool(gene_stats_df):
    """Build the candidate pool with log-transformed matching features.

    Returns (candidate_pool DataFrame, s1 length table, pool-level means/stds).
    """
    gf_geom = pd.read_csv(DATA / "gene_embedding_geometry.csv")
    gf_geom = gf_geom[~gf_geom["gene"].isin(["<pad>", "<mask>", "<cls>", "<eos>"])]

    all_genes = gf_geom[["gene", "ensembl_id", "token_id"]].merge(
        gene_stats_df, on="ensembl_id", how="inner", suffixes=("", "_stats")
    )
    if "token_id_stats" in all_genes.columns:
        all_genes.drop(columns=["token_id_stats"], inplace=True)
    all_genes["gene_class"] = all_genes["gene"].map(assign_gene_class)

    s1 = table_s1[["gene_symbol", "gene_length_bp"]].dropna()
    s1 = s1.rename(columns={"gene_symbol": "gene",
                             "gene_length_bp": "gene_length"})
    all_genes = all_genes.merge(s1, on="gene", how="left")

    # EXCLUDE rather than impute. Genes with no annotation row cannot be
    # matched on length, and median-imputing them would let the matcher treat
    # an unknown value as an average one, which silently biases the null
    # reported an SMD of 0.070 when the true figure was 0.334. A gene missing
    # a matching covariate is simply not an eligible control.
    _n_before = len(all_genes)
    _n_missing = int(all_genes["gene_length"].isna().sum())
    all_genes = all_genes[all_genes["gene_length"].notna()].copy()
    print(f"    Candidate pool: {_n_before} genes -> {len(all_genes)} "
          f"eligible ({_n_missing} excluded: no gene length in Table S1)")
    if len(all_genes) < 5000:
        raise RuntimeError(
            f"Only {len(all_genes)} eligible candidates remain; the pool is "
            f"too small to draw matched control sets.")

    # Every treatment gene must itself be eligible; otherwise the treatment
    # and the control pool are drawn from different universes.
    _tr = pd.read_csv(OUT / "E2_treatment_genes.csv")["gene"]
    _lost = sorted(set(_tr) - set(all_genes["gene"]))
    if _lost:
        raise RuntimeError(
            f"{len(_lost)} treatment gene(s) lack a gene length and were "
            f"excluded from the eligible pool: {_lost[:10]}. Treatment and "
            f"controls would not be comparable.")

    # Log-transform expression and length for matching (Fix 4)
    all_genes["log_expr"] = np.log1p(all_genes["expr_mean"])
    all_genes["log_length"] = np.log1p(all_genes["gene_length"])

    # Compute pool-level z-score parameters (fixed across all draws)
    match_cols = ["log_expr", "expr_breadth", "log_length"]
    pool_stats = {}
    for col in match_cols:
        pool_stats[col] = {
            "mean": float(all_genes[col].mean()),
            "std": float(all_genes[col].std()),
        }

    return all_genes, s1, pool_stats


def _get_treatment_properties(treatment_df, gene_stats_df, s1, pool_stats):
    """Pre-compute log-transformed, z-scored treatment gene properties."""
    treat_props = []
    for _, trow in treatment_df.iterrows():
        match = gene_stats_df.loc[
            gene_stats_df["ensembl_id"] == trow["ensembl_id"]]
        raw_expr = float(match["expr_mean"].values[0]) if len(match) > 0 \
            else 0.0
        raw_breadth = float(match["expr_breadth"].values[0]) if len(match) > 0 \
            else 0.5
        length_vals = s1.loc[s1["gene"] == trow["gene"], "gene_length"].values
        raw_length = float(length_vals[0]) if len(length_vals) > 0 else \
            np.exp(pool_stats["log_length"]["mean"]) - 1

        treat_props.append({
            "gene": trow["gene"],
            "ensembl_id": trow["ensembl_id"],
            "token_id": int(trow["token_id"]),
            "gene_class": trow["gene_class"],
            "expr_mean": raw_expr,
            "expr_breadth": raw_breadth,
            "gene_length": raw_length,
            # z-scored matching features
            "z_log_expr": (np.log1p(raw_expr) - pool_stats["log_expr"]["mean"])
                          / max(pool_stats["log_expr"]["std"], 1e-10),
            "z_breadth": (raw_breadth - pool_stats["expr_breadth"]["mean"])
                         / max(pool_stats["expr_breadth"]["std"], 1e-10),
            "z_log_length": (np.log1p(raw_length)
                             - pool_stats["log_length"]["mean"])
                            / max(pool_stats["log_length"]["std"], 1e-10),
        })
    return treat_props


@functools.lru_cache(maxsize=1)
def _token_class_map():
    """token_id -> mutually exclusive gene class over the GF vocabulary."""
    geom = pd.read_csv(DATA / "gene_embedding_geometry.csv")
    geom = geom[~geom["gene"].isin(["<pad>", "<mask>", "<cls>", "<eos>"])]
    return {int(t): assign_gene_class(g)
            for t, g in zip(geom["token_id"], geom["gene"])}


def _treatment_list_sha256(treatment_df):
    """Stable identity of a treatment set: sha256 of its sorted Ensembl IDs.

    Recorded in every control cache's sidecar and re-derived on reuse, so a
    null drawn against a different treatment set cannot be picked up by a
    later run under the same filename.
    """
    joined = "\n".join(sorted(treatment_df["ensembl_id"]))
    return hashlib.sha256(joined.encode()).hexdigest()


@functools.lru_cache(maxsize=2)
def _expected_draw_profile(sensitivity_arm):
    """Size, class-count vector and treatment tokens a draw must match."""
    treat = get_treatment_genes()
    if sensitivity_arm:
        treat = treat[~treat["gene_class"].isin(
            ["ribosomal", "mitochondrial"])]
    return (len(treat),
            collections.Counter(treat["gene_class"]),
            frozenset(treat["token_id"].astype(int)))


def _controls_fingerprint(cache_path):
    """sha256[:16] of a matched-control cache, used to bind checkpoints."""
    return hashlib.sha256(pathlib.Path(cache_path).read_bytes()).hexdigest()[:16]


def _assert_control_spec(cache_path):
    """Refuse a matched-control cache built under a different specification.

    Control caches have stable filenames, so a cache drawn under a different
    ribosomal panel, or at a different draw count, would
    otherwise be reused silently and the null would not correspond to the
    specification in force. Every field written by build_matched_controls is
    checked, not just the pattern.
    """
    cache_path = pathlib.Path(cache_path)
    spec_path = cache_path.with_name(cache_path.stem + "_spec.json")
    hint = (f"Archive and remove the control caches and ablation checkpoints, "
            f"then re-run --setup:\n"
            f"  rm -f {CACHE}/E2_matched_controls_*.json "
            f"{CACHE}/E2_matched_controls_*_spec.json "
            f"{CACHE}/E2_ablation_ckpt_*.json\n"
            f"(Keep E2_baseline_*, E2_*_tokenized.json, E2_*_gene_stats.csv.)")
    if not spec_path.exists():
        raise RuntimeError(
            f"{cache_path.name} has no _spec.json sidecar, so it predates the "
            f"current ribosomal panel and was drawn under a different "
            f"panel.\n{hint}")
    spec = json.loads(spec_path.read_text())
    sensitivity_arm = cache_path.stem.endswith("_no_ribo_mito")
    treat = get_treatment_genes()
    if sensitivity_arm:
        treat = treat[~treat["gene_class"].isin(
            ["ribosomal", "mitochondrial"])]
    dataset_name = (cache_path.stem
                    .replace("E2_matched_controls_", "")
                    .replace("_no_ribo_mito", ""))
    expected = {"ribosomal_panel_sha256": RIBOSOMAL_PANEL_SHA256,
                "n_bootstrap": N_BOOTSTRAP,
                "dataset": dataset_name,
                "arm": "sensitivity" if sensitivity_arm else "primary",
                "n_draws": N_BOOTSTRAP,
                "genes_per_draw": len(treat),
                "treatment_gene_list_sha256":
                    _treatment_list_sha256(treat),
                "control_cache_sha256":
                    hashlib.sha256(cache_path.read_bytes()).hexdigest()}
    bad = {k: (spec.get(k), v) for k, v in expected.items() if spec.get(k) != v}
    if bad:
        detail = "; ".join(f"{k}: cache={got!r} run={want!r}"
                           for k, (got, want) in bad.items())
        raise RuntimeError(
            f"{cache_path.name} specification mismatch -- {detail}.\n"
            f"The cache does not correspond to the treatment set, arm, "
            f"dataset, eligible pool or null size in force, or it has "
            f"been modified since it was written.\n{hint}")
    draws = json.loads(cache_path.read_text())
    if len(draws) != N_BOOTSTRAP:
        raise RuntimeError(
            f"{cache_path.name} holds {len(draws)} draws but N_BOOTSTRAP="
            f"{N_BOOTSTRAP}.\n{hint}")

    # Every draw is checked structurally rather than trusted. The
    # manuscript states that controls are matched exactly on mutually
    # exclusive gene class and drawn without replacement from genes
    # outside the treatment set, so a cache that violates any of those
    # properties does not correspond to the reported null.
    want_n, want_classes, treat_tokens = _expected_draw_profile(
        cache_path.stem.endswith("_no_ribo_mito"))
    tok_class = _token_class_map()
    for i, draw in enumerate(draws):
        if len(draw) != want_n:
            raise RuntimeError(
                f"{cache_path.name} draw {i} holds {len(draw)} genes, "
                f"expected {want_n}.\n{hint}")
        if len(set(draw)) != len(draw):
            raise RuntimeError(
                f"{cache_path.name} draw {i} repeats a gene.\n{hint}")
        shared = treat_tokens.intersection(int(t) for t in draw)
        if shared:
            raise RuntimeError(
                f"{cache_path.name} draw {i} contains {len(shared)} "
                f"treatment gene(s).\n{hint}")
        got = collections.Counter(
            tok_class.get(int(t), "unmapped") for t in draw)
        if got != want_classes:
            raise RuntimeError(
                f"{cache_path.name} draw {i} class composition "
                f"{dict(got)} does not match the treatment "
                f"{dict(want_classes)}.\n{hint}")


def build_matched_controls(treatment_df, gene_stats_df, dataset_name,
                           label=None):
    """Build 200 class-stratified matched control gene sets.

    Matching uses log1p(expr_mean), expr_breadth, log1p(gene_length)
    with pool-level z-scoring (stable distances across draws).
    Each matched-control draw is independent; within a draw, no gene is reused.
    Saves control gene IDs and a balance table for diagnostics.

    Parameters
    ----------
    label : str, optional
        Suffix for saved files (e.g. "no_ribo_mito"). Defaults to None,
        which uses dataset_name alone.
    """
    all_genes, s1, pool_stats = _prepare_candidate_pool(gene_stats_df)

    treatment_ensembl = set(treatment_df["ensembl_id"])
    treat_props = _get_treatment_properties(
        treatment_df, gene_stats_df, s1, pool_stats
    )

    # Z-score the candidate pool features
    candidate_pool = all_genes[
        (~all_genes["ensembl_id"].isin(treatment_ensembl)) &
        (all_genes["expr_mean"].notna())
    ].copy()
    candidate_pool["z_log_expr"] = (
        (candidate_pool["log_expr"] - pool_stats["log_expr"]["mean"])
        / max(pool_stats["log_expr"]["std"], 1e-10)
    )
    candidate_pool["z_breadth"] = (
        (candidate_pool["expr_breadth"] - pool_stats["expr_breadth"]["mean"])
        / max(pool_stats["expr_breadth"]["std"], 1e-10)
    )
    candidate_pool["z_log_length"] = (
        (candidate_pool["log_length"] - pool_stats["log_length"]["mean"])
        / max(pool_stats["log_length"]["std"], 1e-10)
    )

    suffix = f"_{label}" if label else ""
    fname_base = f"E2_matched_controls_{dataset_name}{suffix}"

    np.random.seed(RANDOM_SEED)
    controls_all = []       # list of token-ID lists
    controls_detail = []    # list of lists of {gene, ensembl_id, token_id, ...}

    for boot in range(N_BOOTSTRAP):
        control_set = []
        control_info = []
        used_in_draw = set()

        for ti, tprops in enumerate(treat_props):
            tc = tprops["gene_class"]
            candidates = candidate_pool[
                (candidate_pool["gene_class"] == tc) &
                (~candidate_pool["token_id"].isin(used_in_draw))
            ].copy()

            # Class matching is exact, as reported in the manuscript. An
            # earlier version dropped the class restriction when fewer
            # than three same-class candidates remained, which would have
            # produced a null that silently did not match on class.
            # nsmallest(5) is well defined for one to four candidates, so
            # a thin stratum is not a reason to relax the constraint.
            if candidates.empty:
                raise RuntimeError(
                    f"No unused same-class control candidate remains for "
                    f"class {tc!r} at treatment position {ti}. The "
                    f"eligible pool is too small for exact class "
                    f"matching at this treatment size.")

            # Euclidean distance in z-scored space
            candidates["dist"] = np.sqrt(
                (candidates["z_log_expr"] - tprops["z_log_expr"])**2 +
                (candidates["z_breadth"] - tprops["z_breadth"])**2 +
                (candidates["z_log_length"] - tprops["z_log_length"])**2
            )

            top5 = candidates.nsmallest(5, "dist")
            chosen = top5.sample(1, random_state=boot * 1000 + ti)
            row = chosen.iloc[0]
            chosen_token = int(row["token_id"])
            control_set.append(chosen_token)
            used_in_draw.add(chosen_token)
            control_info.append({
                "gene": row["gene"],
                "ensembl_id": row["ensembl_id"],
                "token_id": chosen_token,
                "gene_class": row["gene_class"],
                "expr_mean": float(row["expr_mean"]),
                "expr_breadth": float(row["expr_breadth"]),
                "gene_length": float(row["gene_length"]),
                "match_dist": float(row["dist"]),
            })

        controls_all.append(control_set)
        controls_detail.append(control_info)

        if (boot + 1) % 50 == 0:
            print(f"    Bootstrap {boot+1}/{N_BOOTSTRAP}")

    # Save token-ID lists (used by ablation loop)
    with open(CACHE / f"{fname_base}.json", "w") as f:
        json.dump(controls_all, f)

    # Save detailed control gene info (for diagnostics)
    with open(CACHE / f"{fname_base}_detail.json", "w") as f:
        json.dump(controls_detail, f)

    # Spec sidecar: the full specification these controls were drawn under.
    # Every field here is re-derived and compared by _assert_control_spec on
    # reuse, so a cache built under a different ribosomal panel, treatment
    # set, arm, dataset or draw count cannot be silently picked up by a
    # later run, and a cache edited after the fact fails its own hash.
    _cache_bytes = (CACHE / f"{fname_base}.json").read_bytes()
    with open(CACHE / f"{fname_base}_spec.json", "w") as f:
        json.dump({"ribosomal_panel_sha256": RIBOSOMAL_PANEL_SHA256,
                   "n_bootstrap": N_BOOTSTRAP,
                   **panel_provenance(),
                   "dataset": dataset_name,
                   "arm": "sensitivity" if label else "primary",
                   "n_draws": len(controls_all),
                   "genes_per_draw": len(treat_props),
                   "treatment_gene_list_sha256":
                       _treatment_list_sha256(treatment_df),
                   "control_cache_sha256":
                       hashlib.sha256(_cache_bytes).hexdigest()},
                  f, indent=2)

    # ── Balance table (Fix 3) ────────────────────────────────────────────
    _write_balance_table(treat_props, controls_detail, fname_base)
    print(f"    Saved {len(controls_all)} matched control sets "
          f"+ balance table")


def _write_balance_table(treat_props, controls_detail, fname_base):
    """Write a balance table comparing treatment vs controls."""
    # Treatment summary
    t_df = pd.DataFrame(treat_props)

    # Aggregate controls across all 200 draws (mean of means)
    flat_controls = [g for draw in controls_detail for g in draw]
    c_df = pd.DataFrame(flat_controls)

    rows = []
    for col, label in [("expr_mean", "expr_mean"),
                        ("expr_breadth", "expr_breadth"),
                        ("gene_length", "gene_length")]:
        rows.append({
            "variable": label,
            "treatment_mean": float(t_df[col].mean()),
            "treatment_sd": float(t_df[col].std()),
            "control_mean": float(c_df[col].mean()),
            "control_sd": float(c_df[col].std()),
            "smd": float(
                (t_df[col].mean() - c_df[col].mean())
                / max(np.sqrt((t_df[col].std()**2 + c_df[col].std()**2) / 2),
                      1e-10)
            ),
        })

    # Class composition
    t_class = t_df["gene_class"].value_counts(normalize=True)
    c_class = c_df["gene_class"].value_counts(normalize=True)
    for gc in sorted(set(t_class.index) | set(c_class.index)):
        rows.append({
            "variable": f"pct_{gc}",
            "treatment_mean": float(t_class.get(gc, 0)),
            "treatment_sd": 0.0,
            "control_mean": float(c_class.get(gc, 0)),
            "control_sd": 0.0,
            "smd": float(t_class.get(gc, 0) - c_class.get(gc, 0)),
        })

    bal = pd.DataFrame(rows)
    bal.to_csv(OUT / f"{fname_base}_balance.csv", index=False)
    print(f"    Balance table:")
    for _, r in bal.iterrows():
        print(f"      {r['variable']:20s}  treat={r['treatment_mean']:.4f}"
              f"  ctrl={r['control_mean']:.4f}  SMD={r['smd']:+.3f}")


# ── Step 2: Baseline ────────────────────────────────────────────────────────

def select_device(allow_mps=None):
    """Pick compute device: CUDA > MPS (Apple Silicon) > CPU.

    MPS is OPT-IN via E2_USE_MPS=1 because Geneformer/BERT on MPS can differ
    numerically from CPU. Validate with --validate-device before trusting a
    long run; acceptance threshold 1e-5, the same bar that cleared the
    length-aware batching rework.
    """
    import torch
    if torch.cuda.is_available():
        return torch.device("cuda")
    if allow_mps is None:
        allow_mps = os.environ.get("E2_USE_MPS", "0") == "1"
    if allow_mps and torch.backends.mps.is_available():
        if not torch.backends.mps.is_built():
            print("  MPS requested but torch not built with MPS; using CPU")
            return torch.device("cpu")
        return torch.device("mps")
    return torch.device("cpu")


def validate_device(n_cells=64, tol=1e-5):
    """Embed the same cells on CPU and MPS; report max absolute difference.

    Run BEFORE committing to a long MPS run. Returns non-zero if the
    difference exceeds tol, so it can gate a shell pipeline.
    """
    import torch
    from transformers import BertModel

    print("=" * 70)
    print("  DEVICE VALIDATION - CPU vs MPS")
    print("=" * 70)

    if not torch.backends.mps.is_available():
        print("  MPS not available on this machine. Nothing to validate.")
        return 0

    tok_path = CACHE / "E2_pbmc3k_tokenized.json"
    if not tok_path.exists():
        print(f"  Need {tok_path.name}; run --setup first.")
        return 1
    with open(tok_path) as f:
        data = json.load(f)
    tokenized = [np.array(t, dtype=np.int64)
                 for t in data["tokenized"][:n_cells]]

    token_dict, _ = load_geneformer_dicts()
    pad_id = token_dict.get("<pad>", 0)

    model = BertModel.from_pretrained(str(MODEL_DIR))
    model.eval()

    print(f"  Embedding {len(tokenized)} cells on CPU...")
    cpu_embs = extract_cls_embeddings(
        model.to(torch.device("cpu")), tokenized, pad_id,
        torch.device("cpu"), FORWARD_BATCH_SIZE)

    print(f"  Embedding {len(tokenized)} cells on MPS...")
    mps_embs = extract_cls_embeddings(
        model.to(torch.device("mps")), tokenized, pad_id,
        torch.device("mps"), FORWARD_BATCH_SIZE)

    diff = np.abs(cpu_embs.astype(np.float64) - mps_embs.astype(np.float64))
    max_d, mean_d = float(diff.max()), float(diff.mean())
    print(f"\n  max |CPU - MPS|  = {max_d:.3e}")
    print(f"  mean |CPU - MPS| = {mean_d:.3e}")
    print(f"  tolerance        = {tol:.1e}")

    num = (cpu_embs * mps_embs).sum(axis=1)
    den = (np.linalg.norm(cpu_embs, axis=1) * np.linalg.norm(mps_embs, axis=1))
    cos = num / np.maximum(den, 1e-12)
    print(f"  min per-cell cosine(CPU, MPS) = {cos.min():.8f}")

    if max_d <= tol:
        print("\n  PASS - MPS embeddings match CPU within tolerance.")
        print("  Safe to run with E2_USE_MPS=1")
        return 0
    print("\n  FAIL - difference exceeds tolerance. Do NOT use MPS for the")
    print("  production run: results would not be comparable to the PBMC3k")
    print("  baseline computed on CPU.")
    return 1


def run_baseline(datasets=None):
    import torch
    from transformers import BertModel

    print("=" * 70)
    print("  STEP 2: Baseline cell embeddings + linear probe")
    print("=" * 70)

    device = select_device()
    print(f"  Device: {device}")

    print(f"  Loading Geneformer V2-104M...")
    model = BertModel.from_pretrained(str(MODEL_DIR))
    model.eval()
    model = model.to(device)

    token_dict, _ = load_geneformer_dicts()
    pad_id = token_dict.get("<pad>", 0)

    for dataset_name, cache_file in [
        ("pbmc3k", "E2_pbmc3k_tokenized.json"),
        ("tabula_sapiens", "E2_tabula_sapiens_tokenized.json"),
    ]:
        if datasets and dataset_name not in datasets:
            continue
        cache_path = CACHE / cache_file
        if not cache_path.exists():
            print(f"\n  Skipping {dataset_name} (not tokenized)")
            continue

        print(f"\n  {dataset_name}:")
        with open(cache_path) as f:
            data = json.load(f)

        tokenized = [np.array(t, dtype=np.int64) for t in data["tokenized"]]
        cell_types = data["cell_types"]

        # Get CLS embeddings. Persistent chunk dir => a crash costs one
        # 2000-cell chunk, not the whole pass. Safe to reuse here because
        # the baseline embeds the SAME unablated matrix every time; delete
        # the dir if the tokenized input ever changes.
        embs = extract_cls_embeddings(
            model, tokenized, pad_id, device, FORWARD_BATCH_SIZE,
            chunk_path=str(CACHE / f"E2_baseline_chunks_{dataset_name}_"
                                   f"{device.type}")
        )
        print(f"    Embeddings: {embs.shape}")

        # Save baseline embeddings (needed for fixed-probe metric)
        np.save(CACHE / f"E2_baseline_embs_{dataset_name}.npy", embs)

        # Evaluate
        metrics = evaluate_probe(embs, cell_types)
        print(f"    Retrained F1:  {metrics['retrained_f1']:.4f}")
        print(f"    Probe ARI:     {metrics['probe_ari']:.4f}")
        print(f"    Cluster ARI:   {metrics['cluster_ari']:.4f}")
        print(f"    Cluster NMI:   {metrics['cluster_nmi']:.4f}")

        result = {
            "dataset": dataset_name,
            "n_cells": len(cell_types),
            "n_types": len(set(cell_types)),
            **{f"baseline_{k}": v for k, v in metrics.items()},
            "environment": environment_fingerprint(device),
            "table_s1": _table_s1_provenance(),
        }
        with open(OUT / f"E2_baseline_{dataset_name}.json", "w") as f:
            json.dump(result, f, indent=2)
        with open(OUT / f"E2_baseline_{dataset_name}_{device.type}.json",
                  "w") as f:
            json.dump(result, f, indent=2)
        print(f"    Device: {device.type}  ->  also saved "
              f"E2_baseline_{dataset_name}_{device.type}.json")

    del model
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


def make_length_aware_batches(tokenized, max_batch, sq_budget):
    """Group cell indices into batches with bounded attention cost.

    Peak memory for a padded batch is driven by its LONGEST member, since
    attention is O(batch x seq_len^2). Batching by cell count alone spikes
    wherever long cells cluster in the input ordering. Sorting by length and
    capping batch x max_len^2 keeps peak memory flat AND removes wasted
    padding, so it is typically faster despite smaller batches on the tail.

    Returns a list of lists of ORIGINAL indices (order restored by caller).
    """
    order = sorted(range(len(tokenized)), key=lambda i: len(tokenized[i]))
    batches, cur, cur_max = [], [], 0
    for i in order:
        L = len(tokenized[i])
        new_max = max(cur_max, L)
        if cur and (len(cur) + 1 > max_batch
                    or (len(cur) + 1) * new_max * new_max > sq_budget):
            batches.append(cur)
            cur, cur_max = [i], L
        else:
            cur.append(i)
            cur_max = new_max
    if cur:
        batches.append(cur)
    return batches


def extract_cls_embeddings(model, tokenized, pad_id, device, batch_size,
                           chunk_path=None):
    """Extract CLS-token embeddings for all cells, returned in INPUT order.

    Memory safety: batches are length-aware and capped by
    ATTENTION_SQ_BUDGET, so a block of very long cells cannot spike into an
    OOM kill (macOS jetsam SIGKILLs silently — no Python traceback).

    Resume: if chunk_path is given, completed chunks persist on disk and are
    skipped on restart, so a crash costs one chunk rather than the whole run.
    Pass chunk_path=None (default) for ablation runs, where each call embeds
    a DIFFERENT ablated matrix and reusing a fixed chunk dir would silently
    serve the previous control's embeddings.
    """
    import gc
    import tempfile
    import torch

    CHUNK_CELLS = 2000  # cells per on-disk chunk

    n_in = len(tokenized)
    if n_in == 0:
        raise RuntimeError("No sequences to embed")

    # Empty sequences were previously dropped silently, which returned fewer
    # rows than input cells and MISALIGNED embeddings against cell_types —
    # a silent, science-breaking bug. Fail loudly instead.
    n_empty = sum(1 for t in tokenized if len(t) == 0)
    if n_empty:
        raise RuntimeError(
            f"{n_empty} empty sequence(s) found. Filter `tokenized` and "
            f"`cell_types` together upstream — dropping them here would "
            f"misalign embeddings against labels."
        )

    batches = make_length_aware_batches(
        tokenized, max_batch=batch_size, sq_budget=ATTENTION_SQ_BUDGET)

    # Group batches into on-disk chunks
    groups, cur, cur_n = [], [], 0
    for b in batches:
        cur.append(b)
        cur_n += len(b)
        if cur_n >= CHUNK_CELLS:
            groups.append(cur)
            cur, cur_n = [], 0
    if cur:
        groups.append(cur)

    is_temp = chunk_path is None
    if is_temp:
        chunk_path = tempfile.mkdtemp(prefix="E2_emb_chunks_")
    else:
        os.makedirs(chunk_path, exist_ok=True)
        # Persistent chunk dirs are keyed only by dataset name, so a run at a
        # DIFFERENT cell count (e.g. after --subsample) would resume from
        # chunks whose stored indices refer to the old, larger matrix and
        # scatter them into a smaller array -> IndexError, or silently wrong
        # rows. Invalidate on any n_cells mismatch.
        manifest = os.path.join(chunk_path, "manifest.json")
        stale = False
        if os.path.exists(manifest):
            try:
                with open(manifest) as mf:
                    _m = json.load(mf)
                prev = _m.get("n_cells")
                prev_dev = _m.get("device")
                stale = (prev != n_in)
                if stale:
                    print(f"      Chunk cache built for {prev} cells but this "
                          f"run has {n_in}; discarding stale chunks.")
            except Exception:
                stale = True
        elif any(f.startswith("chunk_") for f in os.listdir(chunk_path)):
            print("      Chunk cache has no manifest (written by an older "
                  "version); discarding to be safe.")
            stale = True
        if stale:
            for f in os.listdir(chunk_path):
                if f.startswith("chunk_") or f == "manifest.json":
                    os.remove(os.path.join(chunk_path, f))
        with open(manifest, "w") as mf:
            json.dump({"n_cells": n_in, "device": device.type}, mf)

    max_len_all = max(len(t) for t in tokenized)
    print(f"      {n_in} cells -> {len(batches)} batches / "
          f"{len(groups)} chunks (max_len={max_len_all})")

    done = 0
    for gi, group in enumerate(groups):
        emb_f = os.path.join(chunk_path, f"chunk_{gi:04d}.npy")
        idx_f = os.path.join(chunk_path, f"chunk_{gi:04d}_idx.npy")
        n_group = sum(len(b) for b in group)

        if os.path.exists(emb_f) and os.path.exists(idx_f):
            done += n_group
            print(f"      Chunk {gi+1}/{len(groups)} on disk — resuming "
                  f"({done}/{n_in} cells)")
            continue

        g_embs, g_idx = [], []
        with torch.inference_mode():
            for b in group:
                batch = [tokenized[i] for i in b]
                max_len = max(len(t) for t in batch)
                input_ids = np.full((len(batch), max_len), pad_id,
                                    dtype=np.int64)
                for i, t in enumerate(batch):
                    input_ids[i, :len(t)] = t
                # Mask derived from content, not sequence length — ablated
                # sequences contain trailing pads not to be attended to
                attention_mask = (input_ids != pad_id).astype(np.int64)

                input_ids_t = torch.tensor(input_ids, device=device)
                attention_mask_t = torch.tensor(attention_mask, device=device)

                out = model(input_ids=input_ids_t,
                            attention_mask=attention_mask_t)
                # CLS is position 0 — .copy() breaks shared storage with the
                # full (batch, seq_len, hidden) tensor
                cls_embs = (
                    out.last_hidden_state[:, 0, :]
                    .detach().cpu().contiguous().numpy().copy()
                )
                g_embs.append(cls_embs)
                g_idx.extend(b)

                del out, input_ids_t, attention_mask_t

                # MPS does not release its cache the way CUDA does; without

                # an explicit empty_cache the allocator grows monotonically

                # and OOMs on the long-sequence batches at the tail of the

                # length-sorted order.

                if device.type == 'mps':

                    torch.mps.empty_cache()

        chunk_arr = np.concatenate(g_embs, axis=0)
        np.save(emb_f, chunk_arr)
        np.save(idx_f, np.asarray(g_idx, dtype=np.int64))
        done += n_group
        print(f"      Chunk {gi+1}/{len(groups)} saved ({done}/{n_in} cells)")
        del g_embs, chunk_arr
        gc.collect()

    # Reassemble, scattering each chunk back to its ORIGINAL row positions
    result = None
    for gi in range(len(groups)):
        arr = np.load(os.path.join(chunk_path, f"chunk_{gi:04d}.npy"))
        idx = np.load(os.path.join(chunk_path, f"chunk_{gi:04d}_idx.npy"))
        if result is None:
            result = np.zeros((n_in, arr.shape[1]), dtype=arr.dtype)
        if idx.size and (idx.max() >= n_in or idx.min() < 0):
            raise RuntimeError(
                f"Chunk {gi} holds row index {int(idx.max())} but this run has "
                f"only {n_in} cells. The chunk cache at {chunk_path} belongs "
                f"to a different run. Delete that directory and re-run.")
        result[idx] = arr

    if is_temp:
        for gi in range(len(groups)):
            os.remove(os.path.join(chunk_path, f"chunk_{gi:04d}.npy"))
            os.remove(os.path.join(chunk_path, f"chunk_{gi:04d}_idx.npy"))
        os.rmdir(chunk_path)

    if result is None or result.shape[0] != n_in:
        raise RuntimeError(
            f"Row count mismatch: got "
            f"{None if result is None else result.shape[0]}, expected {n_in}")
    return result


def evaluate_probe(embs, cell_types):
    """5-fold CV linear probe: returns dict of metrics.

    Metrics returned:
      retrained_f1   — macro-F1 from retrained probe (stringent test)
      probe_ari      — ARI between true labels and probe predictions
                       (renamed to clarify this is supervised, not clustering)
      cluster_ari    — ARI between true labels and k-means clusters
      cluster_nmi    — NMI between true labels and k-means clusters
    """
    le = LabelEncoder()
    y = le.fit_transform(cell_types)
    n_classes = len(le.classes_)

    # ── Retrained linear probe (5-fold CV) ───────────────────────────────
    skf = StratifiedKFold(n_splits=N_CV_FOLDS, shuffle=True,
                          random_state=RANDOM_SEED)
    f1_scores = []
    all_preds = np.zeros_like(y)

    for train_idx, test_idx in skf.split(embs, y):
        clf = LogisticRegression(max_iter=2000, random_state=RANDOM_SEED,
                                 class_weight="balanced")
        clf.fit(embs[train_idx], y[train_idx])
        preds = clf.predict(embs[test_idx])
        all_preds[test_idx] = preds
        f1_scores.append(f1_score(y[test_idx], preds, average="macro"))

    retrained_f1 = float(np.mean(f1_scores))
    probe_ari = float(adjusted_rand_score(y, all_preds))

    # ── Unsupervised clustering metrics ──────────────────────────────────
    km = KMeans(n_clusters=n_classes, n_init=10, random_state=RANDOM_SEED)
    cluster_labels = km.fit_predict(embs)
    cluster_ari = float(adjusted_rand_score(y, cluster_labels))
    cluster_nmi = float(normalized_mutual_info_score(y, cluster_labels))

    return {
        "retrained_f1": retrained_f1,
        "probe_ari": probe_ari,
        "cluster_ari": cluster_ari,
        "cluster_nmi": cluster_nmi,
    }


def evaluate_fixed_probe(baseline_embs, ablated_embs, cell_types):
    """Fixed-probe metric: train on baseline, predict on ablated.

    Tests whether the representation moves enough to break an existing
    downstream model (no recalibration allowed).

    NOTE: trains on all baseline cells and predicts on the same cells
    after ablation. This is a representation-shift diagnostic, not a
    generalisation metric. The retrained 5-fold F1 is the gate metric.
    """
    le = LabelEncoder()
    y = le.fit_transform(cell_types)

    # Train on ALL baseline embeddings
    clf = LogisticRegression(max_iter=2000, random_state=RANDOM_SEED,
                             class_weight="balanced")
    clf.fit(baseline_embs, y)

    # Predict on ablated embeddings with the frozen classifier
    preds = clf.predict(ablated_embs)
    fixed_f1 = float(f1_score(y, preds, average="macro"))
    fixed_ari = float(adjusted_rand_score(y, preds))
    return {"fixed_f1": fixed_f1, "fixed_probe_ari": fixed_ari}


# ── Step 3: Ablation ────────────────────────────────────────────────────────

def run_ablation(datasets=None, primary_only=False):
    import torch
    from transformers import BertModel

    print("=" * 70)
    print("  STEP 3: Treatment + matched-control ablations")
    print("=" * 70)

    device = select_device()
    print(f"  Device: {device}")

    model = BertModel.from_pretrained(str(MODEL_DIR))
    model.eval()
    model = model.to(device)

    token_dict, _ = load_geneformer_dicts()
    pad_id = token_dict.get("<pad>", 0)

    # Load treatment gene token IDs
    treatment = pd.read_csv(OUT / "E2_treatment_genes.csv")
    treatment_tokens = treatment["token_id"].tolist()
    print(f"  Treatment (full): {len(treatment_tokens)} genes")

    treatment_no_rm = treatment[~treatment["gene_class"].isin(
        ["ribosomal", "mitochondrial"])]
    treatment_no_rm_tokens = treatment_no_rm["token_id"].tolist()
    print(f"  Treatment (no ribo/mito): {len(treatment_no_rm_tokens)} genes")

    for dataset_name in (datasets or ["pbmc3k", "tabula_sapiens"]):
        tok_path = CACHE / f"E2_{dataset_name}_tokenized.json"
        if not tok_path.exists():
            print(f"\n  Skipping {dataset_name} (not tokenized)")
            continue

        print(f"\n  ── {dataset_name} ──")
        with open(tok_path) as f:
            data = json.load(f)
        tokenized = [np.array(t, dtype=np.int64) for t in data["tokenized"]]
        cell_types = data["cell_types"]

        # Load baseline
        with open(OUT / f"E2_baseline_{dataset_name}.json") as f:
            baseline = json.load(f)
        baseline_f1 = baseline["baseline_retrained_f1"]
        print(f"    Baseline F1: {baseline_f1:.4f}")

        # Load baseline embeddings for fixed-probe
        baseline_embs = np.load(
            CACHE / f"E2_baseline_embs_{dataset_name}.npy")

        # ── Treatment ablation (full k=50) ───────────────────────────────
        print(f"    Treatment ablation (k={TOP_K})...")
        treat_metrics = run_single_ablation(
            model, tokenized, cell_types, treatment_tokens,
            pad_id, device, baseline_embs
        )
        print(f"    Treatment retrained F1: "
              f"{treat_metrics['retrained_f1']:.4f} "
              f"(Δ={treat_metrics['retrained_f1'] - baseline_f1:+.4f})")
        print(f"    Treatment fixed F1:     "
              f"{treat_metrics['fixed_f1']:.4f}")

        # ── Sensitivity: no ribo/mito ────────────────────────────────────
        sens_metrics, k_results = None, {}
        if primary_only:
            print("    [--primary-only] skipping sensitivity arm, k-sweep "
                  "and the sensitivity null band")
        else:
            print(f"    Sensitivity ablation (no ribo/mito)...")
            sens_metrics = run_single_ablation(
                model, tokenized, cell_types, treatment_no_rm_tokens,
                pad_id, device, baseline_embs
            )
            print(f"    Sensitivity retrained F1: "
                  f"{sens_metrics['retrained_f1']:.4f} "
                  f"(Δ={sens_metrics['retrained_f1'] - baseline_f1:+.4f})")

            # ── Sensitivity: k=25 and k=100 (reported, not gated) ───────
            for alt_k in [25, 100]:
                alt_treatment = get_treatment_genes(alt_k)
                alt_tokens = alt_treatment["token_id"].tolist()
                alt_metrics = run_single_ablation(
                    model, tokenized, cell_types, alt_tokens,
                    pad_id, device, baseline_embs
                )
                k_results[alt_k] = alt_metrics
                print(f"    k={alt_k}: retrained F1="
                      f"{alt_metrics['retrained_f1']:.4f} "
                      f"(Δ={alt_metrics['retrained_f1'] - baseline_f1:+.4f})")

        # ── Null band: full treatment controls ───────────────────────────
        ctrl_path = CACHE / f"E2_matched_controls_{dataset_name}.json"
        if not ctrl_path.exists():
            print(f"    No matched controls, skipping null band")
            continue

        ctrl_results = _run_control_null(
            model, tokenized, cell_types, ctrl_path,
            pad_id, device, baseline_embs, baseline_f1,
            dataset_name, "full"
        )

        # ── Null band: sensitivity controls (own matched null) ───────────
        sens_ctrl_path = (CACHE /
            f"E2_matched_controls_{dataset_name}_no_ribo_mito.json")
        sens_ctrl_results = None
        if primary_only:
            pass
        elif sens_ctrl_path.exists():
            sens_ctrl_results = _run_control_null(
                model, tokenized, cell_types, sens_ctrl_path,
                pad_id, device, baseline_embs, baseline_f1,
                dataset_name, "no_ribo_mito"
            )
        else:
            print(f"    No sensitivity controls — skipping sensitivity null")

        # ── Save full results ────────────────────────────────────────────
        results = {
            "dataset": dataset_name,
            "baseline": baseline,
            "treatment": treat_metrics,
            "treatment_delta_f1": float(
                treat_metrics["retrained_f1"] - baseline_f1),
            "sensitivity": sens_metrics,
            "sensitivity_delta_f1": (
                None if sens_metrics is None
                else float(sens_metrics["retrained_f1"] - baseline_f1)),
            "primary_only": bool(primary_only),
            "k_sensitivity": {str(k): v for k, v in k_results.items()},
            "control_results_full": ctrl_results,
            "control_results_no_ribo_mito": sens_ctrl_results,
            "environment": environment_fingerprint(device),
            "table_s1": _table_s1_provenance(),
            "n_bootstrap": N_BOOTSTRAP,
        }
        with open(OUT / f"E2_ablation_{dataset_name}.json", "w") as f:
            json.dump(results, f, indent=2)
        print(f"    Saved ablation results")

    del model
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


def _run_control_null(model, tokenized, cell_types, ctrl_path,
                      pad_id, device, baseline_embs, baseline_f1,
                      dataset_name, null_label):
    """Run the control null band for a given set of matched controls."""
    _assert_control_spec(ctrl_path)
    with open(ctrl_path) as f:
        controls = json.load(f)

    print(f"    Running {null_label} null band "
          f"({len(controls)} controls)...")

    ckpt_path = CACHE / f"E2_ablation_ckpt_{dataset_name}_{null_label}.json"
    ctrl_fp = _controls_fingerprint(ctrl_path)
    ctrl_results = []
    start_idx = 0
    if ckpt_path.exists():
        _ck = json.loads(ckpt_path.read_text())
        if not isinstance(_ck, dict) or "controls_sha256" not in _ck:
            raise RuntimeError(
                f"{ckpt_path.name} is in an older checkpoint format and is not "
                f"bound to a control set, so resuming could mix results from "
                f"different nulls. Delete it and restart this band.")
        if _ck["controls_sha256"] != ctrl_fp:
            raise RuntimeError(
                f"{ckpt_path.name} was written against controls "
                f"{_ck['controls_sha256']} but {ctrl_path.name} is now "
                f"{ctrl_fp}. Resuming would splice two different nulls. "
                f"Delete the checkpoint and restart this band.")
        ctrl_results = _ck["results"]
        start_idx = len(ctrl_results)
        print(f"      Resuming from control {start_idx}/{len(controls)} "
              f"(controls {ctrl_fp})")

    import time
    _t0 = time.time()
    _n_done_this_run = 0

    for ci in range(start_idx, len(controls)):
        ctrl_tokens = controls[ci]
        ctrl_metrics = run_single_ablation(
            model, tokenized, cell_types, ctrl_tokens,
            pad_id, device, baseline_embs
        )
        ctrl_results.append({
            "control_idx": ci,
            **ctrl_metrics,
            "delta_f1": float(
                ctrl_metrics["retrained_f1"] - baseline_f1),
        })

        _n_done_this_run += 1

        # Timing: first control is printed immediately so throughput is known
        # within minutes rather than after the first 10-control checkpoint.
        if _n_done_this_run == 1 or (ci + 1) % 10 == 0:
            with open(ckpt_path, "w") as f:
                json.dump({"controls_sha256": ctrl_fp,
                           "results": ctrl_results}, f)
            _per = (time.time() - _t0) / _n_done_this_run
            _left = len(controls) - (ci + 1)
            print(f"      Control {ci+1}/{len(controls)}: "
                  f"F1={ctrl_metrics['retrained_f1']:.4f} "
                  f"(Δ={ctrl_metrics['retrained_f1']-baseline_f1:+.4f}) "
                  f"[{_per/60:.1f} min/control, "
                  f"ETA {_per*_left/3600:.1f} h for this band]")

    # Clean up checkpoint
    if ckpt_path.exists():
        ckpt_path.unlink()

    return ctrl_results


def run_single_ablation(model, tokenized, cell_types, tokens_to_delete,
                        pad_id, device, baseline_embs=None):
    """Ablate genes, extract embeddings, evaluate both probe types.

    Returns dict with retrained_f1, probe_ari, cluster_ari, cluster_nmi,
    fixed_f1, fixed_probe_ari.
    """
    ablated = [delete_gene_tokens(t, tokens_to_delete, pad_id)
               for t in tokenized]

    embs = extract_cls_embeddings(model, ablated, pad_id, device,
                                 FORWARD_BATCH_SIZE)

    # Retrained probe + clustering metrics
    metrics = evaluate_probe(embs, cell_types)

    # Fixed probe (train on baseline, test on ablated)
    if baseline_embs is not None:
        fixed = evaluate_fixed_probe(baseline_embs, embs, cell_types)
        metrics.update(fixed)

    return metrics


# ── Step 4: Evaluate + verdict ───────────────────────────────────────────────

def evaluate(datasets=None):
    """Gate assessment.

    Deltas are recomputed here from absolute retrained_f1 against the PINNED
    standalone baseline, not read from the ablation JSON. The ablation stores
    deltas against whatever baseline was in memory when it ran, and the probe
    baseline moves by ~6e-4 between execution environments, so those stored
    deltas are not safe to quote.

    The inference statistics -- gate z, tail count and empirical p -- compare
    treatment and control F1 directly and are invariant to the baseline. They
    are final regardless of which baseline is pinned. The descriptive deltas
    and null-band coordinates are not, and are labelled accordingly.
    """
    print("=" * 70)
    print("  STEP 4: Gate assessment")
    print("=" * 70)

    # Honour --datasets. Without this an evaluation aimed at one dataset
    # silently rewrites the other's verdict file, which is how a canonical
    # result gets clobbered by an unrelated run.
    for dataset_name in (datasets or ["pbmc3k", "tabula_sapiens"]):
        results_path = OUT / f"E2_ablation_{dataset_name}.json"
        if not results_path.exists():
            print(f"\n  Skipping {dataset_name} (no results)")
            continue

        with open(results_path) as f:
            results = json.load(f)

        run_baseline = results["baseline"]["baseline_retrained_f1"]

        # ── Pinned baseline is the source of truth for descriptive deltas ──
        pinned_path = OUT / f"E2_baseline_{dataset_name}.json"
        baseline_source, pinned_env = "ablation-run (NO STANDALONE FILE)", None
        pinned_baseline = run_baseline
        if pinned_path.exists():
            with open(pinned_path) as f:
                pin = json.load(f)
            pinned_baseline = pin["baseline_retrained_f1"]
            pinned_env = pin.get("environment")
            baseline_source = str(pinned_path.name)

        baseline_agrees = abs(pinned_baseline - run_baseline) < 1e-9

        print(f"\n  ── {dataset_name} ──")
        print(f"    Pinned baseline F1:     {pinned_baseline:.7f}  "
              f"({baseline_source})")
        if not baseline_agrees:
            print(f"    Ablation-run baseline:  {run_baseline:.7f}  "
                  f"(differs by {pinned_baseline - run_baseline:+.2e})")
            print(f"    -> deltas below are rebased onto the pinned value; "
                  f"z / rank / p are unaffected.")
        if pinned_env is None:
            print(f"    WARNING: pinned baseline carries no environment "
                  f"fingerprint. Regenerate it before quoting deltas.")

        def _arm(ctrl_key, treat_key, delta_key):
            ctrl = results.get(ctrl_key, [])
            treat = results.get(treat_key)
            if not ctrl or treat is None:
                return None
            c = np.array([r["retrained_f1"] for r in ctrl], dtype=float)
            t = float(treat["retrained_f1"])
            m = len(c)
            # Inference: baseline-invariant, computed on absolute F1.
            z = float((t - c.mean()) / c.std(ddof=1))
            b = int((c <= t).sum())
            # Add-one empirical p: the controls are a random sample of the
            # possible matched sets, so a Monte-Carlo p must not be able to
            # reach zero (Phipson & Smyth 2010).
            p_addone = (b + 1) / (m + 1)
            # Description: depends on the baseline.
            dc = c - pinned_baseline
            dt = t - pinned_baseline
            return {
                "n_controls": m,
                "treatment_retrained_f1": t,
                "control_mean_retrained_f1": float(c.mean()),
                # Baseline-invariant effect size: the baseline cancels exactly.
                # This is the quantity to lead with -- treatment vs its own
                # matched controls -- rather than either arm's delta against a
                # baseline that moves with the execution device.
                "treatment_minus_control_mean_f1": float(t - c.mean()),
                "treatment_delta_f1": float(dt),
                "control_mean_delta_f1": float(dc.mean()),
                "control_sd_delta_f1": float(dc.std(ddof=1)),
                "null_band_2_5": float(np.percentile(dc, 2.5)),
                "null_band_97_5": float(np.percentile(dc, 97.5)),
                "inside_band": bool(dt >= np.percentile(dc, 2.5)),
                "gate_z": z,
                "n_controls_at_least_as_damaging": b,
                "empirical_p_raw": b / m,
                "empirical_p_addone": p_addone,
            }

        full = _arm("control_results_full", "treatment", "treatment_delta_f1")
        sens = _arm("control_results_no_ribo_mito", "sensitivity",
                    "sensitivity_delta_f1")

        if results.get("primary_only"):
            print("    (run was --primary-only: no sensitivity arm, k-sweep "
                  "or sensitivity null)")

        for label, a in (("FULL treatment", full), ("SENSITIVITY", sens)):
            if a is None:
                continue
            print(f"")
            print(f"    {label}")
            print(f"      treat - ctrl mean   "
                  f"{a['treatment_minus_control_mean_f1']:+.6f}   "
                  f"<- baseline-invariant")
            print(f"      ΔF1 (vs pinned)     {a['treatment_delta_f1']:+.6f}"
                  f"   (baseline-dependent)")
            print(f"      control mean ΔF1    "
                  f"{a['control_mean_delta_f1']:+.6f} "
                  f"(sd {a['control_sd_delta_f1']:.6f})")
            print(f"      95% null band       "
                  f"[{a['null_band_2_5']:+.6f}, {a['null_band_97_5']:+.6f}]")
            print(f"      inside band         {a['inside_band']}")
            print(f"      gate z              {a['gate_z']:+.4f}")
            print(f"      controls >= damage  "
                  f"{a['n_controls_at_least_as_damaging']}/{a['n_controls']}")
            print(f"      empirical p         "
                  f"{a['empirical_p_addone']:.4f} (add-one; "
                  f"raw {a['empirical_p_raw']:.4f})")

        # Gate is the 2.5th-percentile rule of the matched-control band.
        treat_below = full is not None and not full["inside_band"]
        sens_below = sens is not None and not sens["inside_band"]

        treat_fixed = results["treatment"].get("fixed_f1")
        if treat_fixed is not None:
            print(f"")
            print(f"    Treatment fixed-probe F1: {treat_fixed:.4f} "
                  f"(in-sample; not comparable to the CV baseline)")

        treat_clust = results["treatment"].get("cluster_ari")
        bl_clust = results["baseline"].get("baseline_cluster_ari")
        if treat_clust is not None and bl_clust is not None:
            print(f"    Baseline cluster ARI {bl_clust:.4f} -> treatment "
                  f"{treat_clust:.4f} (Δ={treat_clust - bl_clust:+.4f})")

        if treat_below and sens_below:
            gate = "POSITIVE"
            note = ("Treatment ablation degrades cell-type annotation "
                    "significantly more than expression-matched controls "
                    "(below 2.5th percentile of full null). The effect "
                    "survives removal of ribosomal/mitochondrial outliers "
                    "(below 2.5th percentile of sensitivity-specific null). "
                    "Geometry independently predicts task vulnerability.")
        elif treat_below and not sens_below:
            gate = "PARTIAL"
            note = ("Treatment effect is significant vs full null but "
                    "does not survive removal of ribosomal/mitochondrial "
                    "genes when tested against its own matched null. "
                    "The geometric signal is confounded with gene class.")
        elif not treat_below and sens_below:
            gate = "PARTIAL_INVERSE"
            note = ("Full treatment within null, but non-ribo/mito subset "
                    "shows effect. Unexpected — inspect balance tables.")
        else:
            gate = "NULL"
            note = ("The treatment did not show excess disruption relative to "
                    "its matched-control null in this dataset.")
            if results.get("primary_only"):
                note += (" Only the primary arm was run: the "
                         "ribosomal/mitochondrial sensitivity arm was not "
                         "replicated for this dataset.")

        print(f"\n    GATE: {gate}")
        print(f"    {note}")

        verdict = {
            "dataset": dataset_name,
            "gate": gate,
            "note": note,
            "baseline": {
                "pinned_retrained_f1": pinned_baseline,
                "pinned_source": baseline_source,
                "pinned_has_environment": pinned_env is not None,
                "ablation_run_retrained_f1": run_baseline,
                "agrees": baseline_agrees,
                "note": ("Descriptive deltas and null-band coordinates are "
                         "computed against the pinned baseline. gate_z, the "
                         "tail count and empirical p compare treatment and "
                         "control F1 directly and do not depend on it."),
            },
            "full": full,
            "sensitivity": sens,
            "treatment_fixed_f1": treat_fixed,
            "treatment_cluster_ari": treat_clust,
            "baseline_cluster_ari": bl_clust,
            "k_sensitivity": results.get("k_sensitivity", {}),
            "environment_ablation_run": results.get("environment"),
            "environment_pinned_baseline": pinned_env,
            "table_s1": results.get("table_s1"),
            "n_bootstrap": results.get("n_bootstrap"),
        }
        with open(OUT / f"E2_verdict_{dataset_name}.json", "w") as f:
            json.dump(verdict, f, indent=2, default=str)

    print(f"\n  Saved verdict files.")


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="E2: Expression-matched downstream ablation")
    parser.add_argument("--setup", action="store_true")
    parser.add_argument("--baseline", action="store_true")
    parser.add_argument("--ablation", action="store_true")
    parser.add_argument("--evaluate", action="store_true")
    parser.add_argument("--all", action="store_true")
    parser.add_argument(
        "--primary-only", action="store_true",
        help="Run ONLY the top-k treatment and its matched null. Skips the "
             "sensitivity arm, the k-sweep and the sensitivity null band -- "
             "104 of 205 forward-pass sets. Use for an external replication "
             "where only the primary gate is needed.")
    parser.add_argument(
        "--pin-env", action="store_true",
        help="Print the environment fingerprint and exit. Record this "
             "alongside any number destined for a manuscript (H1).")
    parser.add_argument(
        "--validate-device", action="store_true",
        help="Embed sample cells on CPU and MPS; report max |diff| vs 1e-5. "
             "Run this before any long MPS run.")
    parser.add_argument(
        "--subsample", type=int, default=None, metavar="N",
        help="Stratified subsample to N cells per dataset before tokenizing "
             "(preserves cell-type proportions). Used to make Tabula Sapiens "
             "tractable; record the value in the manuscript Methods.")
    parser.add_argument(
        "--n-bootstrap", type=int, default=None, metavar="N",
        help=f"Override the number of matched-control draws "
             f"(default {N_BOOTSTRAP}). "
             "A floor of 100 applies; fewer widens the null band.")
    parser.add_argument(
        "--datasets", nargs="+", default=None,
        metavar="NAME",
        help="Restrict to these datasets (pbmc3k, tabula_sapiens). "
             "Default: both. TS is ~20x the CPU cost of PBMC3k per pass.")
    args = parser.parse_args()

    if args.pin_env:
        print(json.dumps(environment_fingerprint(), indent=2))
        sys.exit(0)

    if args.validate_device:
        sys.exit(validate_device())

    if args.n_bootstrap is not None:
        if args.n_bootstrap < 100:
            print(f"REFUSED: --n-bootstrap {args.n_bootstrap} is below the "
                  f"floor of 100. A widened null band is how a false null "
                    f"is manufactured.")
            sys.exit(2)
        N_BOOTSTRAP = args.n_bootstrap
        print(f"  N_BOOTSTRAP overridden to {N_BOOTSTRAP}")

    if args.subsample is not None:
        SUBSAMPLE_N = args.subsample
        print(f"  Stratified subsample to {SUBSAMPLE_N} cells per dataset")

    if args.all or args.setup:
        setup(args.datasets)
    if args.all or args.baseline:
        run_baseline(args.datasets)
    if args.all or args.ablation:
        run_ablation(args.datasets, primary_only=args.primary_only)
    if args.all or args.evaluate:
        evaluate(args.datasets)

    if not any([args.setup, args.baseline, args.ablation,
                args.evaluate, args.all]):
        parser.print_help()
