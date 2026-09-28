"""
Special-token sensitivity check for the Geneformer geometry screen (camera-ready
v2, 2026-09-16; Methods 2.1).

notebooks/P01_embedding_geometry.ipynb scores the full 20,275-row Geneformer
input-embedding matrix, including the four special tokens (<pad>, <mask>,
<cls>, <eos>), and drops those rows only from the reported outlier set (412
calls, of which <pad> and <eos> are tokens, giving the 410 reported genes).
This script recomputes the four metrics, the z-scores and the anomaly score
with the four rows removed before the centroid, the neighbours and the
standardisation, and compares the two procedures.

Result on 2026-09-15 (numpy 2.x, exact cosine neighbours):
  outliers 410 = 410, symmetric difference 0 genes, top-50 identical,
  Spearman rho of anomaly scores 1.000, largest |dz| on any gene 0.06
  (isolation), i.e. the reported set does not depend on the choice.

Inputs: data/gene_embeddings.npy and data/gene_names.json as written by
notebooks/D01_model_acquisition.ipynb (the raw matrix is not shipped; the
per-gene metrics it yields are data/gene_embedding_geometry.csv). Without the
raw matrix the script only verifies the 412 -> 410 accounting from the CSV.

Usage: python analysis/check_special_tokens.py
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
SPECIAL = {"<pad>", "<mask>", "<cls>", "<eos>"}
K = 10
Z_CUT = 3.0


def geometry(E, k=K, chunk=2000):
    """norm, Euclidean distance to centroid, cosine to centroid, mean cosine
    distance to k nearest neighbours (self excluded); as in P01."""
    c = E.mean(axis=0)
    norm = np.linalg.norm(E, axis=1)
    dist = np.linalg.norm(E - c, axis=1)
    cos = (E @ c) / (np.clip(norm, 1e-10, None) * max(np.linalg.norm(c), 1e-10))
    En = E / np.clip(norm, 1e-10, None)[:, None]
    n = E.shape[0]
    iso = np.empty(n)
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        D = 1.0 - En[s:e] @ En.T
        D[np.arange(e - s), np.arange(s, e)] = np.inf
        iso[s:e] = np.partition(D, k, axis=1)[:, :k].mean(axis=1)
    M = np.stack([norm, dist, cos, iso], axis=1)
    z = (M - M.mean(axis=0)) / M.std(axis=0, ddof=0)
    return z, np.abs(z).max(axis=1)


def spearman(a, b):
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def csv_accounting():
    import csv
    rows = list(csv.DictReader(open(DATA / "gene_embedding_geometry.csv")))
    calls = [r for r in rows if r["is_outlier"] == "True"]
    tok = [r["gene"] for r in calls if r["gene"] in SPECIAL]
    print(f"gene_embedding_geometry.csv: {len(rows)} rows, {len(calls)} calls, "
          f"{len(tok)} of them special tokens {tok}, {len(calls) - len(tok)} reported genes")


def main():
    csv_accounting()
    emb, names = DATA / "gene_embeddings.npy", DATA / "gene_names.json"
    if not (emb.exists() and names.exists()):
        print("raw matrix not present (run notebooks/D01_model_acquisition.ipynb); "
              "recomputation skipped")
        return 0
    X = np.load(emb).astype(np.float64)
    nm = json.load(open(names))
    sp = np.array([n in SPECIAL for n in nm])
    z_all, A_all = geometry(X)                  # as run: tokens scored, then excluded
    out_all = (A_all > Z_CUT) & ~sp
    keep = ~sp
    z_ex, A_ex = geometry(X[keep])              # tokens removed before scoring
    out_ex = np.zeros(len(nm), bool)
    out_ex[np.where(keep)[0]] = A_ex > Z_CUT
    top_all = set(np.argsort(-np.where(sp, -np.inf, A_all))[:50])
    A_ex_full = np.full(len(nm), -np.inf)
    A_ex_full[keep] = A_ex
    top_ex = set(np.argsort(-A_ex_full)[:50])
    print(f"outliers: scored-then-excluded {out_all.sum()}, removed-before-scoring {out_ex.sum()}, "
          f"symmetric difference {(out_all ^ out_ex).sum()}")
    print(f"top-50 overlap {len(top_all & top_ex)}/50; Spearman rho {spearman(A_all[keep], A_ex):.4f}; "
          f"max |dz| per metric {np.abs(z_all[keep] - z_ex).max(axis=0).round(4)}")
    return 0 if (out_all ^ out_ex).sum() == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
