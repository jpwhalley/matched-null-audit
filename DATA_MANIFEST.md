# Data manifest

What is shipped, what must be downloaded, what is generated, and what each
step costs.

## Shipped inputs

Verified against `data/CHECKSUMS.json` by `E1_stage2_esm2.py --verify-shipped`.

| File | SHA-256 (16) | Size | Contents |
|---|---|---|---|
| `data/Table_S1.csv` | `57443b7225229e0b` | 5.9 MB | Per-gene annotation: geometry scores, cross-model outlier-set membership, expression, gene class, constraint and ClinVar |
| `data/gene_embedding_geometry.csv` | `94d814094f5a519c` | 4.4 MB | Geneformer per-gene geometry: four metrics with z-scores |
| `data/scgpt_gene_embedding_geometry.csv` | `2c7ff95b71a5c237` | 12.3 MB | scGPT per-gene geometry |
| `data/sf_gene_embedding_geometry.csv` | `d3a4924c156c8835` | 2.6 MB | scFoundation per-gene geometry |
| `data/ribosomal_panel.csv` | `acd2a7adf03ad9b4` | 7 KB | Pinned ribosomal panel, HGNC groups 728, 729, 646 |
| `data/ribosomal_panel_provenance.json` | `13e36ef7e78f80eb` | 1 KB | Panel source, retrieval date and hash |
| `outputs/E1_esm2_geometry.csv` | `cc3741ad4a65bde4` | 2.7 MB | ESM-2 protein-embedding geometry, same four metrics |

`Table_S1.csv` is the single annotation table behind Tables 1 and 2 and the
class analyses. Its `outlier_class` column records which models call each gene an
outlier — cross-model set membership, not a biological gene class. Scripts
read the outlier sets from it rather than recalling them, so the sets cannot
drift between analyses. Biological gene class is a separate assignment, made
from the ribosomal panel, mitochondrial symbols, constraint and ClinVar.

`data/ribosomal_panel.csv` is loaded through `analysis/_ribosomal_panel.py`,
which refuses to run if the file's hash disagrees with its provenance record.

## Matched-control draws

`cache/E2_matched_controls_<dataset>[_no_ribo_mito].json` hold the control gene
sets behind the deletion nulls, with a `_spec.json` sidecar recording the panel
hash, draw count, dataset, arm, genes per draw and treatment-list hash.
`E2_downstream_ablation.py` refuses a cache whose sidecar disagrees with the
specification in force.

These are shipped so the null bands can be recomputed without redrawing them.
The tokenised cell matrices they index are not shipped; rebuild those with
`E2_downstream_ablation.py --setup`.

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

Model weights are not redistributed here.

## Generated

| File | Produced by |
|---|---|
| `data/gene_embedding_geometry.csv` and the two companions | `P01_embedding_geometry.ipynb` |
| `outputs/E1_esm2_comparison.csv`, `E1_stage2_verdict.json` | `E1_stage2_esm2.py` |
| `outputs/E3_*.{csv,json}` | `E3_outlier_robustness.py` |
| `outputs/E6_*.{csv,json}` | `E6_class_association.py` |
| `outputs/E8_clinvar_adjusted.{csv,json}` | `E8_clinvar_adjusted.py` |
| `outputs/E10_cross_model_agreement.{csv,json}` | `E10_cross_model_agreement.py` |
| `outputs/E2_treatment_genes.csv`, `E2_matched_controls_*_balance.csv` | `E2_downstream_ablation.py --setup` |
| `outputs/E2_baseline_<dataset>.json` | `E2_downstream_ablation.py --baseline` |
| `outputs/E2_ablation_<dataset>.json` | `E2_downstream_ablation.py --ablation` |
| `outputs/E2_verdict_<dataset>.json` | `E2_downstream_ablation.py --evaluate` |
| `outputs/E7_cluster_metric_diagnostic.{csv,json}` | `E7_cluster_metric_diagnostic.py` |
| `outputs/E9_token_occurrence*.{csv,json}` | `E9_token_occurrence_audit.py` |
| `figures/F1_designs.pdf` … `figures/F4_nullband.pdf` | `make_psb_figures.py` |

Floats in `outputs/` are serialised at the fixed precision defined in
`analysis/_precision.py`, so the drift check in `run_all.sh` is not tripped by
library-version noise. Figures suppress PDF timestamps and are byte-identical
across runs on the same rendering stack.

## Reproduction cost

| Stage | Time | Needs |
|---|---|---|
| E3, E6, E8, E10 | seconds to minutes | shipped data only |
| E1 stage 2 (`--verify-shipped`) | seconds | shipped data only |
| E1 stage 1 | ~10 min | network (MyGene) |
| E1 stage 2 (`--all`) | ~2 h | ESM-2 checkpoint, GPU helpful |
| E2 setup | ~5 min | tokenised PBMC3k |
| **E2 ablation** | **~18 h** | Geneformer checkpoint; MPS or CUDA advisable |
| E7, E9, figures | seconds | E2 outputs |
