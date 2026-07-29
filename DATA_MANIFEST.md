# Data manifest

Every input the analyses need, and where it comes from. Three categories:
**shipped** (in this repository), **downloaded** (fetched by the D-series
notebooks from the original source), and **generated** (produced by scripts
here).

Checksums for shipped files are in `data/CHECKSUMS.json` and are also recorded
inside the output JSONs of any script that reads them.

---

## Shipped

| File | SHA-256 (16) | Size | Contents |
|---|---|---|---|
| `data/Table_S1.csv` | `57443b7225229e0b` | 5.9 MB | Per-gene annotation: geometry scores for all three models, pLI, LOEUF, ClinVar status, expression, breadth, gene length, cross-model outlier class. 18,915 genes. |
| `data/gene_embedding_geometry.csv` | `94d814094f5a519c` | 4.4 MB | Geneformer per-gene geometry: four metrics, z-scores, composite anomaly score, outlier flag. 20,275 rows including four special tokens. |
| `outputs/E1_esm2_geometry.csv` | `cc3741ad4a65bde4` | 2.7 MB | ESM-2 protein-embedding geometry, same four metrics. Shipped because regenerating it requires the ESM-2 checkpoint and several GPU-hours. |

### ⚠️ On `Table_S1.csv`

**Use only this copy.** An earlier copy elsewhere in the project carried gene
lengths for 11,752 of 18,915 genes; the rest were median-imputed downstream.
That silently broke the length-matched design. Three numbers tell the story:
the old controls *appeared* balanced at SMD **0.070**; evaluating those same
controls against complete lengths gave **0.334**; rebuilding the controls with
complete lengths restored balance to **0.077**. The same table also moved an
adjusted odds ratio across the significance boundary.

This copy has 18,911/18,915 lengths (99.98%). Scripts assert ≥99% coverage and
refuse to run below it. If you substitute your own annotation table, that
assertion is the thing protecting you.

---

## Downloaded (D-series notebooks; not redistributed)

| Resource | Version / accession | Notebook | Licence |
|---|---|---|---|
| Geneformer | `Geneformer-V2-104M` | `D01_model_acquisition.ipynb` | Apache-2.0, via Hugging Face |
| scGPT | whole-human checkpoint | `D01_model_acquisition.ipynb` | MIT |
| scFoundation | released checkpoint | `D01_model_acquisition.ipynb` | per authors |
| ESM-2 | `esm2_t33_650M_UR50D` | `D01_model_acquisition.ipynb` | MIT |
| MANE Select | v1.3 (pinned) | `D02_external_databases.ipynb` | public domain |
| gnomAD constraint | v4.1 | `D02_external_databases.ipynb` | ODbL |
| ClinVar | as dated in the notebook | `D02_external_databases.ipynb` | public domain |
| PBMC3k | 10x Genomics, via scanpy | `D04_reference_datasets.ipynb` | CC-BY |
| Tabula Sapiens Immune | CellxGene `78b60b70-129a-4a6d-b15f-825b241eec66` | `D04_reference_datasets.ipynb` | CC-BY |

Model weights and single-cell matrices are **not** redistributed here. Pin
versions when re-downloading: MANE Select in particular is versioned, and a
different release changes the gene-to-protein mapping in stage 3.

---

## Generated

| File | Produced by | Notes |
|---|---|---|
| `outputs/E3_*.{csv,json}` | `E3_outlier_robustness.py` | Caller robustness; includes all six callers, three of them degenerate and excluded with reasons. |
| `outputs/E1_stage2_verdict.json` | `E1_stage2_esm2.py` | ESM-2 comparison verdict. |
| `outputs/E6_*.{csv,json}` | `E6_class_association.py` | Class associations under both schemes; per-class scFM-only breakdown. |
| `outputs/E8_clinvar_adjusted.{csv,json}` | `E8_clinvar_adjusted.py` | Covariate-adjusted disease association, full specification grid. |
| `outputs/E2_treatment_genes.csv` | `E2_downstream_ablation.py --setup` | The top-50 treatment set with class assignments. |
| `outputs/E2_ablation_pbmc3k.json` | `E2_downstream_ablation.py --ablation` | 200 draws against the corrected length-matched controls. |
| `outputs/E2_baseline_pbmc3k.json` | `E2_downstream_ablation.py --baseline` | Pinned MPS baseline, device-matched to the ablation; carries an `environment` fingerprint. `_cpu` and `_mps` copies retained alongside. |
| `outputs/E2_verdict_pbmc3k.json` | `E2_downstream_ablation.py --evaluate` | Gate z, tail count and add-one empirical p, recomputed from absolute F1 against the pinned baseline. |
| `outputs/E7_*`, `outputs/E9_*` | `E7_*.py`, `E9_*.py` | Regenerated against the corrected control sets. |
| `figures/F4_nullband_pbmc3k.pdf` | `make_psb_figures.py` | Corrected controls; shows z = −1.56 and z = +0.35, both inside. |
| `cache/` | various | Regenerable intermediates; git-ignored. |

**E2, E7, E9 and Figure 4 form one consistent set and were added together.**
While the corrected rerun was outstanding they were deliberately absent rather
than stale: a repository containing a mixture of generations is worse than one
that is visibly incomplete.

---

## Reproduction cost

| Stage | Time | Needs |
|---|---|---|
| E3, E6, E8 | seconds to minutes | shipped data only |
| E1 stage 1 | ~10 min | network (MyGene) |
| E1 stage 2 | ~2 h | ESM-2 checkpoint, GPU helpful |
| E2 setup | ~5 min | tokenised PBMC3k |
| **E2 ablation** | **~18 h** | Geneformer checkpoint; MPS or CUDA advisable |
| E7, E9, figures | seconds | E2 outputs |

`SKIP_ABLATION=1 ./run_all.sh` reproduces everything except the deletion test.
