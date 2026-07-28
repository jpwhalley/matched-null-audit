"""
E1 Stage 2 — ESM-2 sequence-embedding control.

Generates ESM-2 embeddings for canonical proteins mapped to each scFM's
gene vocabulary, runs the identical four-metric geometric audit, and
compares outlier sets against the scFM results.

Canonical protein assignment: MANE-prioritised gene-to-protein mapping
(MANE Select v1.3, pinned for reproducibility) with
reviewed UniProt (Swiss-Prot) sequence retrieval. One sequence per gene.
Provenance (UniProt accession actually embedded) is recorded per gene.

Usage:
  # Step 1: Build canonical mapping + fetch sequences (can run anywhere)
  python E1_stage2_esm2.py --fetch-sequences

  # Step 2: Run ESM-2 inference (needs GPU or large-RAM machine)
  python E1_stage2_esm2.py --run-esm2

  # Step 3: Audit + compare (runs after embeddings are cached)
  python E1_stage2_esm2.py --audit

  # Or run everything:
  python E1_stage2_esm2.py --all

Requirements:
  pip install torch fair-esm pandas scipy scikit-learn mygene

Outputs (in revision/outputs/ unless noted):
  E1_canonical_proteins.csv — one-protein-per-gene mapping table
  E1_protein_sequences.fasta — canonical sequences (revision/cache/)
  E1_esm2_embeddings.npy — mean-pooled ESM-2 vectors (n_genes × 1280)
  E1_esm2_gene_order.json — gene order matching the embedding matrix
  E1_esm2_geometry.csv — four-metric audit on ESM-2 space
  E1_esm2_comparison.csv — per-model outlier overlap statistics
  E1_stage2_verdict.json — structured verdict
"""

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

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

# Load Table_S1 for gene classes
table_s1 = pd.read_csv(DATA / "Table_S1.csv")


# ── Shared helpers ───────────────────────────────────────────────────────────

def assign_gene_class(sym: str) -> str:
    import re
    sym_u = sym.upper()
    if sym_u.startswith("MT-"):
        return "mitochondrial"
    if re.match(r"^(RPL|RPS|MRPL|MRPS)\d", sym_u):
        return "ribosomal"
    constrained = set(
        table_s1.loc[(table_s1["pLI"] > 0.9) | (table_s1["LOEUF"] < 0.35),
                     "gene_symbol"].str.upper()
    )
    if sym_u in constrained:
        return "constrained"
    disease = set(
        table_s1.loc[table_s1["clinvar_disease"] == True,
                     "gene_symbol"].str.upper()
    )
    if sym_u in disease:
        return "disease"
    return "other"


# ── Step 1: Canonical mapping + sequence fetch ───────────────────────────────

def fetch_sequences():
    import mygene

    print("=" * 70)
    print("  STEP 1: Canonical protein mapping + sequence fetch")
    print("=" * 70)

    # Load MANE Select
    mane_path = CACHE / "MANE.GRCh38.v1.3.summary.txt.gz"
    if not mane_path.exists():
        print("  Downloading MANE Select v1.3...")
        url = ("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/"
               "release_1.3/MANE.GRCh38.v1.3.summary.txt.gz")
        urllib.request.urlretrieve(url, mane_path)

    mane = pd.read_csv(mane_path, sep="\t", compression="gzip")
    mane = mane[mane["MANE_status"] == "MANE Select"].copy()
    mane["ens_base"] = mane["Ensembl_Gene"].str.replace(r"\.\d+$", "",
                                                         regex=True)
    print(f"  MANE Select entries: {len(mane)}")

    # Build lookup tables
    ens_to_refseq = dict(zip(mane["ens_base"], mane["RefSeq_prot"]))
    ens_to_ensp   = dict(zip(mane["ens_base"], mane["Ensembl_prot"]))
    ens_to_symbol = dict(zip(mane["ens_base"], mane["symbol"]))
    sym_to_ens    = dict(zip(mane["symbol"], mane["ens_base"]))
    sym_to_refseq = dict(zip(mane["symbol"], mane["RefSeq_prot"]))

    # Load mygene cache for fallback Swiss-Prot IDs
    mg_cache_path = CACHE / "E1_mygene_cache.json"
    if mg_cache_path.exists():
        with open(mg_cache_path) as f:
            mg_cache = json.load(f)
    else:
        mg_cache = {}

    # Load feasibility table for protein-annotated genes
    feas = pd.read_csv(OUT / "E1_feasibility_table.csv")

    # Build canonical mapping: union of all protein-coding genes across models
    # Use MANE Select RefSeq protein as primary, Swiss-Prot as fallback
    all_genes = set()
    gene_to_ensembl = {}

    for model in ["Geneformer", "scGPT", "scFoundation"]:
        mdf = feas[(feas["model"] == model) & (feas["has_protein"] == True)]
        for _, row in mdf.iterrows():
            sym = row["gene_symbol"]
            ens = row.get("ensembl_gene", None)
            if pd.notna(sym):
                all_genes.add(sym)
            if pd.notna(ens):
                gene_to_ensembl[sym] = ens

    print(f"  Unique protein-annotated genes across models: {len(all_genes)}")

    # Assign canonical protein per gene
    canonical = []
    for sym in sorted(all_genes):
        ens = gene_to_ensembl.get(sym, sym_to_ens.get(sym, None))

        # Try MANE Select first
        refseq_prot = None
        source = None
        if ens and ens in ens_to_refseq:
            refseq_prot = ens_to_refseq[ens]
            source = "MANE_Select"
        elif sym in sym_to_refseq:
            refseq_prot = sym_to_refseq[sym]
            source = "MANE_Select_by_symbol"
            ens = sym_to_ens.get(sym, ens)

        # Fallback to Swiss-Prot from mygene
        swissprot = None
        if refseq_prot is None:
            # Check mygene cache across all model queries
            for model_key in ["geneformer", "scgpt", "scfoundation"]:
                mg_results = mg_cache.get(model_key, {})
                # Try by Ensembl ID or by symbol
                for query_key in [ens, sym]:
                    if query_key and query_key in mg_results:
                        r = mg_results[query_key]
                        up = r.get("uniprot", {})
                        if isinstance(up, dict):
                            sp = up.get("Swiss-Prot", None)
                            if isinstance(sp, list):
                                sp = sp[0]
                            if sp:
                                swissprot = sp
                                source = "UniProt_SwissProt"
                                break
                if swissprot:
                    break

        canonical.append({
            "gene_symbol": sym,
            "ensembl_gene": ens,
            "refseq_protein": refseq_prot,
            "uniprot_id": swissprot,
            "canonical_source": source or "unmapped",
            "protein_id": refseq_prot or swissprot,
            "gene_class": assign_gene_class(sym),
        })

    canon_df = pd.DataFrame(canonical)
    n_mapped = (canon_df["canonical_source"] != "unmapped").sum()
    print(f"  Canonical mapping: {n_mapped}/{len(canon_df)} genes mapped "
          f"({100*n_mapped/len(canon_df):.1f}%)")
    print(f"    MANE Select: "
          f"{(canon_df['canonical_source'].str.startswith('MANE')).sum()}")
    print(f"    Swiss-Prot fallback: "
          f"{(canon_df['canonical_source'] == 'UniProt_SwissProt').sum()}")
    print(f"    Unmapped: "
          f"{(canon_df['canonical_source'] == 'unmapped').sum()}")

    canon_df.to_csv(OUT / "E1_canonical_proteins.csv", index=False)

    # Fetch sequences from UniProt using batched OR queries
    print(f"\n  Fetching protein sequences from UniProt...")

    fasta_path = CACHE / "E1_protein_sequences.fasta"
    seq_cache_path = CACHE / "E1_sequence_cache.json"

    if seq_cache_path.exists():
        with open(seq_cache_path) as f:
            seq_cache = json.load(f)
        print(f"    Loaded {len(seq_cache)} cached sequences")
    else:
        seq_cache = {}

    mapped = canon_df[canon_df["canonical_source"] != "unmapped"]
    to_fetch_syms = [r["gene_symbol"] for _, r in mapped.iterrows()
                     if r["gene_symbol"] not in seq_cache]

    print(f"    Need to fetch: {len(to_fetch_syms)} sequences")

    def _parse_fasta_multi(fasta_text):
        """Parse multi-entry FASTA, return {gene_name: (accession, sequence)}.

        Extracts both GN= (gene name) and the UniProt accession from the
        header so we can record provenance — which exact protein was embedded.
        """
        results = {}
        current_gene = None
        current_acc = None
        current_seq = []
        for line in fasta_text.split("\n"):
            if line.startswith(">"):
                if current_gene and current_seq:
                    results[current_gene] = (
                        current_acc, "".join(current_seq))
                current_seq = []
                current_gene = None
                current_acc = None
                # Header format: >sp|ACCESSION|ENTRY_NAME ...  GN=SYMBOL ...
                parts = line.split()
                if parts:
                    acc_fields = parts[0].lstrip(">").split("|")
                    if len(acc_fields) >= 2:
                        current_acc = acc_fields[1]
                    else:
                        current_acc = acc_fields[0]
                for part in parts:
                    if part.startswith("GN="):
                        current_gene = part[3:]
                        break
            else:
                current_seq.append(line.strip())
        if current_gene and current_seq:
            results[current_gene] = (current_acc, "".join(current_seq))
        return results

    # Batch query: groups of 80 genes per OR query (URL length safe)
    # seq_cache stores {gene: {"acc": accession, "seq": sequence}}
    batch_size = 80
    for i in range(0, len(to_fetch_syms), batch_size):
        batch = to_fetch_syms[i:i+batch_size]
        gene_clauses = " OR ".join(f"gene_exact:{s}" for s in batch)
        query = f"({gene_clauses}) AND organism_id:9606 AND reviewed:true"
        url = (f"https://rest.uniprot.org/uniprotkb/stream?"
               f"query={urllib.parse.quote(query)}&format=fasta")
        try:
            req = urllib.request.Request(
                url,
                headers={"Accept": "text/plain",
                         "User-Agent": "Python/E1_stage2"}
            )
            resp = urllib.request.urlopen(req, timeout=60)
            fasta_text = resp.read().decode()
            parsed = _parse_fasta_multi(fasta_text)
            # Only take genes we asked for (avoid partial symbol matches)
            batch_set = set(batch)
            for gene, (acc, seq) in parsed.items():
                if gene in batch_set and gene not in seq_cache:
                    seq_cache[gene] = {"acc": acc, "seq": seq}
        except Exception as e:
            print(f"    Batch {i//batch_size} failed: {e}")
            # Fall back to individual queries for this batch
            for sym in batch:
                if sym in seq_cache:
                    continue
                try:
                    furl = (f"https://rest.uniprot.org/uniprotkb/search?"
                            f"query=gene_exact:{sym}+AND+organism_id:9606"
                            f"+AND+reviewed:true&format=fasta&size=1")
                    req2 = urllib.request.Request(
                        furl,
                        headers={"Accept": "text/plain",
                                 "User-Agent": "Python/E1_stage2"}
                    )
                    resp2 = urllib.request.urlopen(req2, timeout=10)
                    ft = resp2.read().decode()
                    if ft.strip() and ">" in ft:
                        lines = ft.strip().split("\n")
                        # Extract accession from header
                        hdr = lines[0]
                        acc_parts = hdr.lstrip(">").split("|")
                        acc = acc_parts[1] if len(acc_parts) >= 2 \
                            else acc_parts[0]
                        seq = "".join(l for l in lines[1:]
                                      if not l.startswith(">"))
                        if seq:
                            seq_cache[sym] = {"acc": acc, "seq": seq}
                except Exception:
                    pass

        # Save checkpoint every 10 batches
        if (i // batch_size) % 10 == 0:
            with open(seq_cache_path, "w") as f:
                json.dump(seq_cache, f)
            print(f"    Progress: {min(i+batch_size, len(to_fetch_syms))}"
                  f"/{len(to_fetch_syms)} "
                  f"({len(seq_cache)} cached)")

        # Rate limit: UniProt allows ~100 req/sec for streaming
        time.sleep(0.5)

    # Save sequence cache
    with open(seq_cache_path, "w") as f:
        json.dump(seq_cache, f)

    # Helper: extract sequence string from cache (handles old plain-string
    # format and new {acc, seq} dict format)
    def _get_seq(entry):
        if isinstance(entry, dict):
            return entry["seq"]
        return entry  # old format: plain string

    def _get_acc(entry):
        if isinstance(entry, dict):
            return entry.get("acc", "unknown")
        return "unknown"

    # Write FASTA file with accession in header
    with open(fasta_path, "w") as f:
        for sym in sorted(seq_cache.keys()):
            acc = _get_acc(seq_cache[sym])
            f.write(f">{sym}|{acc}\n{_get_seq(seq_cache[sym])}\n")

    n_with_seq = sum(1 for _, r in canon_df.iterrows()
                     if r["gene_symbol"] in seq_cache)
    print(f"\n  Sequences obtained: {n_with_seq}/{len(mapped)} mapped genes "
          f"({100*n_with_seq/len(mapped):.1f}%)")

    # Sequence length statistics
    lengths = [len(_get_seq(seq_cache[sym])) for sym in seq_cache]
    print(f"  Sequence lengths: median={np.median(lengths):.0f}, "
          f"mean={np.mean(lengths):.0f}, max={max(lengths)}, "
          f"n>1022={sum(1 for l in lengths if l > 1022)}")

    # Add sequence info + provenance to canonical table
    canon_df["has_sequence"] = canon_df["gene_symbol"].isin(seq_cache)
    canon_df["seq_length"] = canon_df["gene_symbol"].map(
        lambda s: len(_get_seq(seq_cache[s])) if s in seq_cache else np.nan
    )
    canon_df["uniprot_embedded"] = canon_df["gene_symbol"].map(
        lambda s: _get_acc(seq_cache[s]) if s in seq_cache else None
    )
    canon_df.to_csv(OUT / "E1_canonical_proteins.csv", index=False)

    print(f"\n  Saved: E1_canonical_proteins.csv, E1_protein_sequences.fasta")


# ── Step 2: ESM-2 inference ─────────────────────────────────────────────────

def run_esm2():
    import torch
    import esm

    print("=" * 70)
    print("  STEP 2: ESM-2 embedding inference")
    print("=" * 70)

    # Load sequences (handles both old plain-string and new {acc,seq} format)
    seq_cache_path = CACHE / "E1_sequence_cache.json"
    with open(seq_cache_path) as f:
        seq_cache_raw = json.load(f)

    def _get_seq(entry):
        return entry["seq"] if isinstance(entry, dict) else entry

    # Load canonical mapping to get only mapped genes
    canon = pd.read_csv(OUT / "E1_canonical_proteins.csv")
    mapped_genes = sorted(
        canon.loc[canon["gene_symbol"].isin(seq_cache_raw), "gene_symbol"]
    )
    print(f"  Genes to embed: {len(mapped_genes)}")

    # Load ESM-2
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"  Device: {device}")
    print(f"  Loading esm2_t33_650M_UR50D...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval().to(device)
    batch_converter = alphabet.get_batch_converter()

    # Resume from checkpoint if available
    ckpt_path = CACHE / "E1_esm2_checkpoint.npz"
    embeddings = {}
    if ckpt_path.exists():
        ckpt = np.load(ckpt_path, allow_pickle=True)
        ckpt_genes = list(ckpt["genes"])
        ckpt_embs = ckpt["embeddings"]
        for g, e in zip(ckpt_genes, ckpt_embs):
            embeddings[g] = e
        print(f"  Resumed from checkpoint: {len(embeddings)} genes done")

    # Embed proteins
    ESM2_CONTEXT = 1022  # max tokens excluding BOS/EOS
    CHECKPOINT_EVERY = 500

    for idx, sym in enumerate(mapped_genes):
        if sym in embeddings:
            continue  # already done (checkpoint resume)

        seq = _get_seq(seq_cache_raw[sym])

        if len(seq) <= ESM2_CONTEXT:
            # Single pass
            data = [(sym, seq)]
            _, _, tokens = batch_converter(data)
            tokens = tokens.to(device)
            with torch.no_grad():
                results = model(tokens, repr_layers=[33])
            # Mean-pool residue representations, excluding BOS (0) and EOS (-1)
            rep = results["representations"][33][0, 1:-1, :].cpu().numpy()
            emb = rep.mean(axis=0)
        else:
            # Windowed mean-pool for long proteins
            window_size = ESM2_CONTEXT
            stride = window_size // 2  # 50% overlap
            window_embs = []
            for start in range(0, len(seq), stride):
                end = min(start + window_size, len(seq))
                subseq = seq[start:end]
                if len(subseq) < 10:
                    continue
                data = [(f"{sym}_{start}", subseq)]
                _, _, tokens = batch_converter(data)
                tokens = tokens.to(device)
                with torch.no_grad():
                    results = model(tokens, repr_layers=[33])
                rep = results["representations"][33][0, 1:-1, :].cpu().numpy()
                window_embs.append(rep.mean(axis=0))
                if start + window_size >= len(seq):
                    break
            emb = np.mean(window_embs, axis=0)

        embeddings[sym] = emb

        if (idx + 1) % 100 == 0 or idx == 0:
            print(f"    {idx+1}/{len(mapped_genes)} "
                  f"({sym}, len={len(seq)})")

        # Checkpoint every CHECKPOINT_EVERY new embeddings
        n_done = len(embeddings)
        if n_done % CHECKPOINT_EVERY == 0:
            g_list = sorted(embeddings.keys())
            np.savez(ckpt_path,
                     genes=np.array(g_list, dtype=object),
                     embeddings=np.array([embeddings[g] for g in g_list]))
            print(f"    [checkpoint saved: {n_done} genes]")

    # Final save
    gene_order = sorted(embeddings.keys())
    emb_matrix = np.array([embeddings[g] for g in gene_order])
    np.save(OUT / "E1_esm2_embeddings.npy", emb_matrix)
    with open(OUT / "E1_esm2_gene_order.json", "w") as f:
        json.dump(gene_order, f)

    # Clean up checkpoint
    if ckpt_path.exists():
        ckpt_path.unlink()

    print(f"\n  Saved: E1_esm2_embeddings.npy ({emb_matrix.shape})")
    print(f"         E1_esm2_gene_order.json ({len(gene_order)} genes)")


# ── Step 3: Audit + compare ──────────────────────────────────────────────────

def audit_and_compare():
    print("=" * 70)
    print("  STEP 3: Geometric audit + scFM comparison")
    print("=" * 70)

    # Load ESM-2 embeddings
    emb = np.load(OUT / "E1_esm2_embeddings.npy")
    with open(OUT / "E1_esm2_gene_order.json") as f:
        gene_order = json.load(f)

    print(f"  ESM-2 embeddings: {emb.shape}")

    # ── Four-metric audit ────────────────────────────────────────────────
    # 1. Norm
    norms = np.linalg.norm(emb, axis=1)
    # 2. Distance from centroid
    centroid = emb.mean(axis=0)
    dists = np.linalg.norm(emb - centroid, axis=1)
    # 3. Cosine similarity to centroid
    cos_sims = np.array([
        np.dot(emb[i], centroid) / (np.linalg.norm(emb[i]) * np.linalg.norm(centroid))
        for i in range(len(emb))
    ])
    # 4. Isolation score (mean cosine distance to k nearest neighbours)
    # Chunked to avoid allocating a full n×n matrix (~2.7 GB for 19k genes)
    from sklearn.metrics.pairwise import cosine_similarity
    k = 10
    n = len(emb)
    isolation_scores = np.zeros(n)
    chunk_size = min(1000, n)  # process 1000 rows at a time
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        # (end-start) × n similarity matrix — manageable
        chunk_sims = cosine_similarity(emb[start:end], emb)
        for local_i in range(end - start):
            global_i = start + local_i
            sims = chunk_sims[local_i]
            # Exclude self
            sims[global_i] = -np.inf
            # k nearest = k highest cosine similarities
            top_k = np.partition(sims, -k)[-k:]
            isolation_scores[global_i] = 1.0 - top_k.mean()

    # Z-scores
    def zscore(x):
        return (x - x.mean()) / x.std()

    norm_z = zscore(norms)
    dist_z = zscore(dists)
    cos_z  = zscore(cos_sims)
    iso_z  = zscore(isolation_scores)

    # Anomaly score = max |z| across metrics
    anomaly_scores = np.maximum.reduce([
        np.abs(norm_z), np.abs(dist_z), np.abs(cos_z), np.abs(iso_z)
    ])

    # Outlier calling: |z| > 3 on any metric (same as scFM method)
    is_outlier = (np.abs(norm_z) > 3) | (np.abs(dist_z) > 3) | \
                 (np.abs(cos_z) > 3) | (np.abs(iso_z) > 3)

    # Build geometry table
    esm2_geom = pd.DataFrame({
        "gene": gene_order,
        "norm": norms,
        "dist_from_centroid": dists,
        "cos_to_centroid": cos_sims,
        "isolation_score": isolation_scores,
        "norm_zscore": norm_z,
        "dist_zscore": dist_z,
        "cos_zscore": cos_z,
        "isolation_zscore": iso_z,
        "anomaly_score": anomaly_scores,
        "is_outlier": is_outlier,
        "gene_class": [assign_gene_class(g) for g in gene_order],
    })
    esm2_geom.to_csv(OUT / "E1_esm2_geometry.csv", index=False)

    compare_against_scfms(esm2_geom)


def verify_shipped():
    """Recompute the scFM comparison from the SHIPPED ESM-2 geometry.

    The geometry table is shipped because regenerating it needs the ESM-2
    checkpoint and several GPU-hours, and the 97 MB embedding array is not
    committed. This path verifies the shipped file against
    data/CHECKSUMS.json and then recomputes everything downstream of it, so
    the comparison and verdict are genuinely regenerated rather than assumed.
    """
    import hashlib

    geom_path = OUT / "E1_esm2_geometry.csv"
    if not geom_path.exists():
        raise FileNotFoundError(
            f"{geom_path} not found. Either restore the shipped file or run "
            f"--all to regenerate it from the ESM-2 checkpoint.")

    checks = REPO / "data" / "CHECKSUMS.json"
    if checks.exists():
        expected = json.load(open(checks)).get(
            "outputs/E1_esm2_geometry.csv", {}).get("sha256_16")
        actual = hashlib.sha256(geom_path.read_bytes()).hexdigest()[:16]
        if expected and actual != expected:
            raise RuntimeError(
                f"Checksum mismatch for {geom_path.name}: expected "
                f"{expected}, got {actual}. The shipped geometry has been "
                f"modified; regenerate with --all or restore the committed "
                f"file.")
        print(f"  Shipped ESM-2 geometry verified: sha256[:16]={actual}")
    else:
        print("  WARNING: data/CHECKSUMS.json absent; cannot verify.")

    esm2_geom = pd.read_csv(geom_path)
    print(f"  Loaded {len(esm2_geom)} rows from shipped geometry")
    compare_against_scfms(esm2_geom)


def compare_against_scfms(esm2_geom):
    """Everything downstream of the ESM-2 geometry table."""
    gene_order = esm2_geom["gene"].tolist()
    esm2_outliers = set(esm2_geom.loc[esm2_geom["is_outlier"], "gene"])
    print(f"  ESM-2 outliers: {len(esm2_outliers)}")

    # ── Compare against each scFM ────────────────────────────────────────
    MODEL_FILES = {
        "Geneformer":   DATA / "gene_embedding_geometry.csv",
        "scGPT":        DATA / "scgpt_gene_embedding_geometry.csv",
        "scFoundation": DATA / "sf_gene_embedding_geometry.csv",
    }

    # Load E3 robust cores
    with open(OUT / "E3_robust_core.json") as f:
        robust_cores = json.load(f)

    comparison_rows = []
    esm2_gene_set = set(gene_order)

    print(f"\n  Comparing ESM-2 vs scFM outliers (shared protein-coding subset):")

    for model_name, geom_path in MODEL_FILES.items():
        scfm = pd.read_csv(geom_path)
        if model_name == "Geneformer":
            scfm = scfm[~scfm["gene"].isin(["<pad>", "<mask>", "<cls>", "<eos>"])]

        sc_col = ("anomaly_score_with_isolation"
                  if "anomaly_score_with_isolation" in scfm.columns
                  else "anomaly_score")

        # Restrict to shared genes (protein-coding, in ESM-2 space)
        shared = set(scfm["gene"]) & esm2_gene_set
        scfm_shared = scfm[scfm["gene"].isin(shared)].copy()
        esm2_shared = esm2_geom[esm2_geom["gene"].isin(shared)].copy()

        scfm_outliers = set(scfm_shared.loc[scfm_shared["is_outlier"], "gene"])
        esm2_out_shared = esm2_outliers & shared

        # Overlap statistics
        intersection = scfm_outliers & esm2_out_shared
        union = scfm_outliers | esm2_out_shared
        jaccard = len(intersection) / len(union) if union else 0

        # Hypergeometric test
        M = len(shared)        # population
        n = len(scfm_outliers) # successes in population
        N = len(esm2_out_shared)  # draws
        k = len(intersection)  # observed successes
        hyper_p = 1.0 - stats.hypergeom.cdf(k - 1, M, n, N)

        # Spearman correlation of anomaly scores
        scfm_for_merge = scfm_shared[["gene", sc_col]].rename(
            columns={sc_col: "anom_scfm"})
        esm2_for_merge = esm2_shared[["gene", "anomaly_score"]].rename(
            columns={"anomaly_score": "anom_esm2"})
        merged = scfm_for_merge.merge(esm2_for_merge, on="gene")
        rho, rho_p = stats.spearmanr(merged["anom_scfm"], merged["anom_esm2"])

        # Category concordance
        scfm_out_classes = scfm_shared.loc[
            scfm_shared["is_outlier"],
            "gene"
        ].map(assign_gene_class).value_counts()
        esm2_out_classes = esm2_shared.loc[
            esm2_shared["is_outlier"],
            "gene"
        ].map(assign_gene_class).value_counts()

        print(f"\n  {model_name} (shared={len(shared)}):")
        print(f"    scFM outliers in shared set: {len(scfm_outliers)}")
        print(f"    ESM-2 outliers in shared set: {len(esm2_out_shared)}")
        print(f"    Intersection: {len(intersection)}")
        print(f"    Jaccard: {jaccard:.3f}")
        print(f"    Hypergeometric p: {hyper_p:.2e}")
        print(f"    Spearman rho: {rho:.3f} (p={rho_p:.2e})")

        print(f"    scFM outlier classes:  "
              f"{dict(scfm_out_classes.head(5))}")
        print(f"    ESM-2 outlier classes: "
              f"{dict(esm2_out_classes.head(5))}")

        if intersection:
            overlap_classes = pd.Series(
                [assign_gene_class(g) for g in intersection]
            ).value_counts()
            print(f"    Overlap classes:       {dict(overlap_classes)}")
            print(f"    Overlap genes:         {sorted(intersection)}")

        # Top-50 overlap (ranked by anomaly score in each space)
        scfm_top50 = set(
            scfm_shared.nlargest(50, sc_col)["gene"]
        )
        esm2_top50 = set(
            esm2_shared.nlargest(50, "anomaly_score")["gene"]
        )
        top50_overlap = scfm_top50 & esm2_top50
        print(f"    Top-50 overlap: {len(top50_overlap)}/50")

        # Category enrichment comparison (Fisher's for each class)
        categories = ["ribosomal", "mitochondrial", "constrained", "disease"]
        cat_concordance = {}
        for cat in categories:
            scfm_cat_in  = sum(1 for g in scfm_outliers
                               if assign_gene_class(g) == cat)
            scfm_cat_out = sum(1 for g in (shared - scfm_outliers)
                               if assign_gene_class(g) == cat)
            esm2_cat_in  = sum(1 for g in esm2_out_shared
                               if assign_gene_class(g) == cat)
            esm2_cat_out = sum(1 for g in (shared - esm2_out_shared)
                               if assign_gene_class(g) == cat)

            _, scfm_p = stats.fisher_exact(
                [[scfm_cat_in,
                  len(scfm_outliers) - scfm_cat_in],
                 [scfm_cat_out,
                  len(shared) - len(scfm_outliers) - scfm_cat_out]],
                alternative="greater"
            )
            _, esm2_p = stats.fisher_exact(
                [[esm2_cat_in,
                  len(esm2_out_shared) - esm2_cat_in],
                 [esm2_cat_out,
                  len(shared) - len(esm2_out_shared) - esm2_cat_out]],
                alternative="greater"
            )
            cat_concordance[cat] = {
                "scfm_enriched": scfm_p < 0.05,
                "esm2_enriched": esm2_p < 0.05,
                "concordant": (scfm_p < 0.05) == (esm2_p < 0.05),
            }

        n_concordant = sum(1 for v in cat_concordance.values()
                           if v["concordant"])
        print(f"    Category concordance: {n_concordant}/{len(categories)}")

        comparison_rows.append({
            "model": model_name,
            "n_shared": len(shared),
            "n_scfm_outliers": len(scfm_outliers),
            "n_esm2_outliers": len(esm2_out_shared),
            "n_intersection": len(intersection),
            "jaccard": jaccard,
            "hypergeom_p": hyper_p,
            "spearman_rho": rho,
            "spearman_p": rho_p,
            "top50_overlap": len(top50_overlap),
            "cat_concordance": n_concordant,
            "overlap_genes": sorted(intersection),
        })

    comp_df = pd.DataFrame(comparison_rows)
    comp_df.to_csv(OUT / "E1_esm2_comparison.csv", index=False)

    # ── Verdict ──────────────────────────────────────────────────────────
    mean_jaccard = comp_df["jaccard"].mean()
    mean_rho = comp_df["spearman_rho"].mean()

    if mean_jaccard < 0.15 and mean_rho < 0.3:
        interpretation = "DIVERGENT"
        note = ("Outlier sets are largely different between scFM and ESM-2 "
                "embeddings. Geometric outliers are model/tokenisation-specific, "
                "not intrinsic biological properties.")
    elif mean_jaccard > 0.5 and mean_rho > 0.6:
        interpretation = "CONVERGENT"
        note = ("Same genes are geometric outliers in both expression-trained "
                "and sequence-trained embeddings. Outlier status reflects "
                "intrinsic biological properties, not model-specific failure.")
    else:
        interpretation = "PARTIAL"
        note = ("Partial overlap: some gene classes (likely sequence-intrinsic) "
                "are shared outliers, while expression-exposure-driven outliers "
                "are scFM-specific. Supports a mixed interpretation.")

    verdict = {
        "interpretation": interpretation,
        "note": note,
        "mean_jaccard": float(mean_jaccard),
        "mean_spearman": float(mean_rho),
        "n_esm2_outliers_total": int(len(esm2_outliers)),
        "per_model": [row.to_dict() for _, row in comp_df.iterrows()],
        "prediction_check": {
            "predicted": "DIVERGENT (different outliers, limited overlap)",
            "predicted_by": "Justin P. Whalley",
            "predicted_date": "2026-07-17",
        },
    }

    with open(OUT / "E1_stage2_verdict.json", "w") as f:
        json.dump(verdict, f, indent=2, default=str)

    print(f"\n\n{'='*70}")
    print(f"  E1 STAGE 2 VERDICT: {interpretation}")
    print(f"  {note}")
    print(f"  Mean Jaccard: {mean_jaccard:.3f}")
    print(f"  Mean Spearman: {mean_rho:.3f}")
    print(f"{'='*70}")

    print(f"\n  Saved:")
    print(f"    E1_esm2_geometry.csv")
    print(f"    E1_esm2_comparison.csv")
    print(f"    E1_stage2_verdict.json")


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="E1 Stage 2: ESM-2 control")
    parser.add_argument("--fetch-sequences", action="store_true",
                        help="Step 1: build canonical mapping + fetch sequences")
    parser.add_argument("--run-esm2", action="store_true",
                        help="Step 2: run ESM-2 inference (needs GPU)")
    parser.add_argument("--audit", action="store_true",
                        help="Step 3: geometric audit + comparison")
    parser.add_argument("--all", action="store_true",
                        help="Run all steps (needs ESM-2 checkpoint + GPU)")
    parser.add_argument("--verify-shipped", action="store_true",
                        help="Verify the shipped ESM-2 geometry against "
                             "data/CHECKSUMS.json and recompute the scFM "
                             "comparison and verdict from it. No checkpoint "
                             "or GPU required.")
    args = parser.parse_args()

    if args.all or args.fetch_sequences:
        fetch_sequences()
    if args.all or args.run_esm2:
        run_esm2()
    if args.verify_shipped:
        verify_shipped()
    if args.all or args.audit:
        audit_and_compare()

    if not any([args.fetch_sequences, args.run_esm2, args.audit, args.all,
                args.verify_shipped]):
        parser.print_help()
