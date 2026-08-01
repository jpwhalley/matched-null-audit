# Traceability

Each reported analysis, the script that produces it, the inputs it reads and
the files it writes. Values are not repeated here; read them from the outputs.

`local` regenerates from shipped inputs alone. `checkpoint` needs model weights
or the ~18 h ablation and is shipped so the reported numbers can be checked
against the artefact.

## Figures

| Manuscript item | Script | Inputs | Outputs | Repro |
|---|---|---|---|---|
| Figure 1, model designs and shared screen | `make_psb_figures.py::fig1_designs` | none (schematic) | `figures/F1_designs.pdf` | local |
| Figure 2, caller robustness | `make_psb_figures.py::fig2_stability` | `outputs/E3_calibrated_summary.csv` | `figures/F2_stability.pdf` | local |
| Figure 3, ESM-2 comparison | `make_psb_figures.py::fig3_esm2` | `outputs/E6_scfm_only_by_class.csv` | `figures/F3_esm2.pdf` | local |
| Figure 4, matched deletion nulls (3 panels) | `make_psb_figures.py::fig4_nullband` | `outputs/E2_ablation_<ds>.json`, `outputs/E2_baseline_<ds>.json` | `figures/F4_nullband.pdf` | local from shipped E2 outputs |

## Tables

| Manuscript item | Script | Inputs | Outputs | Repro |
|---|---|---|---|---|
| Table 1, cross-model agreement | `E10_cross_model_agreement.py` | `data/Table_S1.csv` | `outputs/E10_cross_model_agreement.{csv,json}` | local |
| Table 2a, covariate-adjusted ClinVar | `E8_clinvar_adjusted.py` | `data/Table_S1.csv` | `outputs/E8_clinvar_adjusted.{csv,json}` | local |
| Table 2b, unadjusted class association | `E6_class_association.py` | `data/Table_S1.csv`, `data/gene_embedding_geometry.csv`, `data/ribosomal_panel.csv` | `outputs/E6_class_association.{csv,json}` | local |

## Methods

| Manuscript item | Script | Inputs | Outputs | Repro |
|---|---|---|---|---|
| §2.1 geometry screen, four metrics | `notebooks/P01_embedding_geometry.ipynb` | model checkpoints | `data/*gene_embedding_geometry.csv` | checkpoint |
| §2.3 MANE mapping for ESM-2 | `E1_stage1_mapping.py` | MANE Select v1.3, UniProt | gene-to-protein mapping | checkpoint |
| §2.4 treatment and control construction | `E2_downstream_ablation.py --setup` | `data/Table_S1.csv`, tokenised cells | `outputs/E2_treatment_genes.csv`, `outputs/E2_matched_controls_*_balance.csv`, `cache/E2_matched_controls_*.json` | checkpoint |
| §2.4 pinned baseline and device | `E2_downstream_ablation.py --baseline` | tokenised cells, Geneformer | `outputs/E2_baseline_<ds>.json` | checkpoint |
| §2.5 ribosomal panel definition | `_ribosomal_panel.py` | `data/ribosomal_panel.csv`, `..._provenance.json` | consumed by E1, E2, E3, E6 and E8 | local |
| §2.5 class schemes, overlapping and exclusive | `E6_class_association.py` | `data/Table_S1.csv` | `outputs/E6_class_association.{csv,json}` | local |

## Results

| Manuscript item | Script | Inputs | Outputs | Repro |
|---|---|---|---|---|
| §3.1 outlier counts per model | `E10_cross_model_agreement.py` | `data/Table_S1.csv` | `outputs/E10_cross_model_agreement.json` | local |
| §3.1 pairwise overlap, Jaccard, score correlation | `E10_cross_model_agreement.py` | as above | as above | local |
| §3.2 caller stability and robust cores | `E3_outlier_robustness.py` | the three geometry CSVs, `data/Table_S1.csv` | `outputs/E3_calibrated_summary.csv`, `E3_enrichment_full.csv`, `E3_robust_core.json`, `E3_degenerate_diagnostics.csv`, `E3_gate_verdict.json` | local |
| §3.3 ESM-2 recurrence and per-class breakdown | `E1_stage2_esm2.py --verify-shipped` | `outputs/E1_esm2_geometry.csv`, `data/Table_S1.csv` | `outputs/E1_esm2_comparison.csv`, `E1_stage2_verdict.json` | local |
| §3.3 ESM-2 geometry itself | `E1_stage2_esm2.py --all` | ESM-2 checkpoint | `outputs/E1_esm2_geometry.csv` | checkpoint |
| §3.4 deletion contrast, gate statistic, null band | `E2_downstream_ablation.py --evaluate` | `outputs/E2_ablation_<ds>.json`, `cache/E2_matched_controls_<ds>*.json` | `outputs/E2_verdict_<ds>.json` | checkpoint |
| §3.4 sensitivity arm | as above, `no_ribo_mito` cache and `sensitivity` key | as above | as above | checkpoint |
| §3.4 exclusion of clustering metrics | `E7_cluster_metric_diagnostic.py` | `outputs/E2_ablation_<ds>.json`, `E2_baseline_<ds>.json` | `outputs/E7_cluster_metric_diagnostic.{csv,json}` | local |
| §3.4 token-level exposure | `E9_token_occurrence_audit.py` | `cache/E2_<ds>_tokenized.json`, `cache/E2_matched_controls_<ds>*.json`, `outputs/E2_treatment_genes.csv` | `outputs/E9_token_occurrence*.{csv,json}` | checkpoint |
| §3.4 matching balance, standardised mean differences | `E2_downstream_ablation.py --setup` | as §2.4 | `outputs/E2_matched_controls_*_balance.csv` | checkpoint |
| §3.4 Tabula Sapiens replication | `E2_downstream_ablation.py --datasets tabula_sapiens` | as above | `outputs/E2_{ablation,verdict,baseline}_tabula_sapiens.json` | checkpoint |
| §3.5 adjusted ClinVar association, all models | `E8_clinvar_adjusted.py` | `data/Table_S1.csv` | `outputs/E8_clinvar_adjusted.{csv,json}` | local |
| §3.5 scheme dependence of class enrichment | `E6_class_association.py` | as Table 2b | `outputs/E6_class_association.{csv,json}`, `E6_scfm_only_by_class.csv` | local |

## Environment

Every result file carries an environment fingerprint recording library
versions, thread settings and compute device. `pyproject.toml` and `uv.lock`
pin the environment; `data/CHECKSUMS.json` pins the shipped inputs.
