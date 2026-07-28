# Manuscript–Code Traceability Matrix

> **AI-assisted documentation.** This matrix and the cross-references between
> manuscript claims and code outputs were assembled with the assistance of
> Claude (Anthropic) via Cowork. All scientific analyses, interpretations and
> editorial decisions were made by the author. Readers can use this document to
> verify any claim in the manuscript against the code and saved output that
> produced it.

Every figure panel, table entry and in-text statistic maps to a script, its
saved output file, and the script or notebook that displays it. (Numbered
display notebooks are not yet written; figures are currently built directly by
`analysis/make_psb_figures.py`.)

**Status:** Figures 1–3 and §3.5 are final. Figure 4 and the §3.4 values marked
⏳ are pending the corrected length-matched deletion rerun, and will be filled
from the regenerated `E2_ablation_pbmc3k.json`.

All numbers below are produced from the **shipped** `data/Table_S1.csv`
(SHA-256 `57443b7225229e0b`, 18,911/18,915 gene lengths). Earlier drafts of this
matrix quoted counts from a length-sparse copy of the same table; those are
superseded.

---

## Figure 1 — The five-stage workflow

**Built by:** `analysis/make_psb_figures.py::fig1_workflow`
**Output:** `figures/F1_workflow.pdf`

Schematic only; no data dependencies. Stage outcomes shown on the figure are
summaries of the results below.

---

## Figure 2 — Caller robustness

**Built by:** `analysis/make_psb_figures.py::fig2_stability`
**Data:** `outputs/E3_calibrated_summary.csv` ← `analysis/E3_outlier_robustness.py`

| Claim | Value | Source column |
|---|---|---|
| Geneformer containment under MAD z>3 | 410/410 = 1.000 | `containment` |
| Geneformer Spearman ρ | 0.956 | `spearman_rho` |
| Geneformer top-50 overlap | 50/50 | `top50_overlap` |
| Geneformer Jaccard (set-size artefact) | 0.492 | `jaccard` |
| scFoundation containment | 159/164 = 0.970 | `containment` |
| scFoundation ρ | 0.999 | `spearman_rho` |
| scGPT containment | 156/188 = 0.830 | `containment` |
| scGPT top-50 overlap | 22/50 | `top50_overlap` |
| scFoundation MAD z>3.5 containment saturated | 27/164 = 0.165 | `n_outliers`, `containment` |
| Geneformer IQR containment saturated | 43/410 = 0.105 | `n_outliers`, `containment` |
| Degenerate: GMM minority component | 24.20–49.17% | `outputs/E3_degenerate_diagnostics.csv` |
| Degenerate: percentile 1/99 union rates | 4.96 / 6.77 / 5.92% | `outputs/E3_jaccard_results.csv` |

---

## Figure 3 — Non-transcriptomic (ESM-2) control

**Built by:** `analysis/make_psb_figures.py::fig3_esm2`
**Data:** `outputs/E6_scfm_only_by_class.csv` ← `analysis/E6_class_association.py`
**Upstream:** `analysis/E1_stage2_esm2.py` → `outputs/E1_esm2_geometry.csv`

| Claim | Value | Source |
|---|---|---|
| Mean Jaccard, scFM vs ESM-2 | 0.023 | `outputs/E1_stage2_verdict.json` |
| Mean Spearman | 0.042 | `outputs/E1_stage2_verdict.json` |
| Geneformer ∩ ESM-2 | 49 of 389/570 | `E1_stage2_verdict.json::per_model` |
| Geneformer top-50 overlap | 4 | `E1_stage2_verdict.json::top50_overlap` |
| Shared gene universe | 19,017 | `E6_scfm_only_by_class.csv::n_shared_genes` |
| Ribosomal scFM-only | 47/68 | `E6_scfm_only_by_class.csv` |
| Mitochondrial scFM-only | 2/10 | `E6_scfm_only_by_class.csv` |
| Constrained scFM-only | 86/87 | `E6_scfm_only_by_class.csv` |
| Disease scFM-only | 79/82 | `E6_scfm_only_by_class.csv` |

**Two restrictions are required** to reproduce the constrained figure: the
19,017-gene shared universe, and mutually-exclusive class precedence
(mitochondrial → ribosomal → constrained → disease → other). Omitting the second
gives 139 rather than 87.

---

## Figure 4 — The matched deletion test (claim gate)

**Built by:** `analysis/make_psb_figures.py::fig4_nullband`
**Data:** `outputs/E2_ablation_pbmc3k.json` ← `analysis/E2_downstream_ablation.py`

| Claim | Value | Source |
|---|---|---|
| Baseline macro-F1 (pinned) | 0.9247927 | `outputs/E2_baseline_pbmc3k.json` |
| Treatment retrained F1 | ⏳ | `E2_ablation_pbmc3k.json::treatment.retrained_f1` |
| Treatment ΔF1 | ⏳ | derived |
| Control Δ mean, 95% band | ⏳ | derived from `control_results_full` |
| Gate statistic z | ⏳ | derived |
| One-sided empirical p | ⏳ | derived |
| Controls at least as damaging | ⏳ | derived |
| Sensitivity (36 genes) ΔF1, z | ⏳ | `sensitivity`, `control_results_no_ribo_mito` |
| k-sweep ΔF1 (k=25/50/100) | ⏳ | `k_sensitivity` |
| Treatment composition | 19 constrained · 13 ribosomal · 9 other · 8 disease · 1 mitochondrial | `outputs/E2_treatment_genes.csv` |

Matching balance — `outputs/E2_matched_controls_pbmc3k_balance.csv`:

| Variable | Treatment | Control | SMD |
|---|---|---|---|
| expr_mean | 6.4263 | 3.9514 | **0.288** |
| expr_breadth | 0.4215 | 0.4123 | 0.019 |
| gene_length | 528,680 | 474,122 | 0.077 |
| class proportions | — | — | 0.000 (exact) |

Expression shows **limited common support**: only 28 of 18,283 eligible
candidates reach the treatment's 75th expression percentile, and 10 reach the
90th, against a treatment maximum of 46.1 versus a pool 99.9th percentile of
13.0. This is limited common support: the imbalance reflects scarcity of comparable candidates, not matcher failure.

Excluded from the headline figure — `analysis/E7_cluster_metric_diagnostic.py`:

| Metric | % of controls improving | 95% band as % of baseline | Usable |
|---|---|---|---|
| retrained macro-F1 | 45% | 0.7% | yes |
| cluster NMI | 60% | 10.7% | no |
| cluster ARI | 64% | 42.7% | no |

Token-level exposure — `analysis/E9_token_occurrence_audit.py`:

| Arm | Treatment | Controls | z |
|---|---|---|---|
| Full (50 genes) | 21.07 tokens/cell | 20.62 | +1.73 |
| Sensitivity (36 genes) | 9.42 | 9.29 | +0.62 |

Token occurrence does not indicate treatment underexposure.

---

## Table 1 / §3.5 — Covariate-aware disease association

**Data:** `outputs/E8_clinvar_adjusted.{csv,json}` ← `analysis/E8_clinvar_adjusted.py`

Universe: 18,911 complete cases (388 Geneformer outliers, 72 shared).
Non-mitochondrial sensitivity: 18,898 (378, 64).

| Tier | Gene set | Unadjusted OR | Adjusted OR | 95% CI | p |
|---|---|---|---|---|---|
| **PRIMARY** Firth, all genes | all Geneformer | 1.16 | **1.05** | 0.82–1.35 | 0.680 |
| **PRIMARY** Firth, all genes | shared (exploratory) | 3.59 | 1.45 | 0.80–2.71 | 0.222 |
| SENSITIVITY logistic, non-mito | all Geneformer | 1.10 | 1.05 | 0.82–1.35 | 0.682 |
| SENSITIVITY logistic, non-mito | shared (exploratory) | 3.05 | 1.47 | 0.80–2.71 | 0.215 |
| DIAGNOSTIC mito covariate omitted | shared | 3.59 | 2.28 | 1.26–4.13 | 0.007 |

The DIAGNOSTIC row is **mis-specified** and reported only to document what
happens when a perfectly separated class (all 13 MT- genes are ClinVar-positive)
is left unadjusted. Its estimates must not be quoted.

Class-scheme comparison — `outputs/E6_class_association.csv`:

| Gene set | Scheme | Disease OR | p |
|---|---|---|---|
| shared GF∩scGPT | overlapping | 3.594 | 5.0e-07 |
| shared GF∩scGPT | mutually exclusive | 0.385 | 1.6e-03 |
| all GF outliers | overlapping | 1.181 | 0.104 |
| all GF outliers | mutually exclusive | 0.538 | 1.7e-07 |

The two schemes estimate **different quantities**, not the same quantity two
ways: the mutually-exclusive `disease` class excludes genes already assigned to
constrained, ribosomal or mitochondrial. 1,784 of 8,287 ClinVar genes (22%) are
also constrained.

Geneformer class associations (original |z|>3 call) —
`outputs/E3_enrichment_full.csv`:

| Class | OR | Direction |
|---|---|---|
| mitochondrial | 165.48 | enriched |
| ribosomal | 36.71 | enriched |
| constrained | 1.59 | enriched |
| disease (residual class) | 0.54 | negative — **depends on mutually exclusive class precedence** |

---

## Methods statistics

| Claim | Value | Source |
|---|---|---|
| Geneformer vocabulary | 20,271 genes | `E3_degenerate_diagnostics.csv::total_n` |
| scGPT vocabulary | 60,694 | same |
| scFoundation vocabulary | 19,264 | same |
| Geneformer outliers | 410 | `E3_calibrated_summary.csv` |
| Annotation table coverage | 18,911/18,915 lengths | `E8_clinvar_adjusted.json::provenance` |
| Annotation table SHA-256 | `57443b7225229e0b` | same |
| Canonical pinned environment | Python 3.11.15, numpy 1.26.4, sklearn 1.8.0, threads=1 | `E2_baseline_pbmc3k.json::environment` |

**On the 412 → 410 discrepancy.** The preprint reported 412 Geneformer outliers
across 20,275 genes; this work reports 410 across 20,271. The difference is the
exclusion of four special tokens (`<pad>`, `<mask>`, `<cls>`, `<eos>`), two of
which were geometric outliers. Special tokens are not genes and are excluded
from every analysis here.
