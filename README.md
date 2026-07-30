# Comparative geometry of gene embeddings in single-cell foundation models

Code and reproducibility materials for:

> **Comparative geometry of gene embeddings in single-cell foundation models:
> model-specific outliers and limited predictive value in Geneformer**
> Whalley, J.P. (2026). Submitted to *Pacific Symposium on Biocomputing 2027*.

Preprint: [doi:10.64898/2026.06.22.733850](https://doi.org/10.64898/2026.06.22.733850)

---

## What this repository is for

Every figure, table and in-text statistic in the manuscript maps to a specific
script or notebook cell here, and to the saved output it came from. `MANUSCRIPT_TRACEABILITY.md` is that mapping. A
reader who wants to check a number should be able to find it in one lookup
without running anything.

Notebooks in `notebooks/` are committed **with their outputs executed**, so the
analysis can be read without a GPU, without the model checkpoints, and without
re-running the expensive steps.

> **Status.** The D-series and `P01` are executed and current. The numbered
> display notebooks (one per figure) are **not yet written** — figures are built
> directly by `analysis/make_psb_figures.py` in the meantime. All analysis
> outputs, including the matched deletion test, are present and current.

---

## Layout

```
prespecification/   dated analysis plan, fixed before the deletion compute
analysis/           canonical implementations (the scripts that produced results)
notebooks/          executed notebooks: D-series acquire, P-series analyse,
                    numbered series display results and build figures
figures/            manuscript figures as PDF
outputs/            saved result files (CSV/JSON) consumed by the notebooks
```

The current release pipeline runs left to right:

```
D-series / P01  ->  analysis/E*.py  ->  outputs/*.{csv,json}  ->  analysis/make_psb_figures.py  ->  figures/
(acquire/screen)    (compute)            (results)                  (display)
```

---

## The five audit stages

The same screen is applied to Geneformer, scGPT and scFoundation, and the
resulting outlier sets are then tested for what they support. Each stage can
falsify a claim independently.

| Stage | Question | Implementation |
|---|---|---|
| 1. Geometry screen | Which genes are geometrically extreme? | `notebooks/P01_embedding_geometry.ipynb` |
| 2. Caller robustness | Is the outlier set an artefact of thresholding? | `analysis/E3_outlier_robustness.py` |
| 3. Non-transcriptomic control | Are the same genes outliers in protein-sequence space? | `analysis/E1_stage1_mapping.py`, `analysis/E1_stage2_esm2.py` |
| 4. Matched deletion test | Does removing them affect a downstream task more than matched controls? | `analysis/E2_downstream_ablation.py` |
| 5. Covariate-aware annotation | Does outlier status carry independent disease signal? Run for all three models under one specification. | `analysis/E8_clinvar_adjusted.py`, `analysis/E6_class_association.py` |

Supporting diagnostics:

| Script | Purpose |
|---|---|
| `E7_cluster_metric_diagnostic.py` | Quantifies outcome precision and documents why macro-F1, rather than clustering metrics, was retained as the headline measure. |
| `E9_token_occurrence_audit.py` | Measures token-level exposure of treatment vs control gene sets. The treatment genes are not under-represented, which rules out one route to the null but not others. |

---

## Pre-specification

`prespecification/analysis_plan_2026-07-07.md` predates all deletion compute. It
fixes class-stratified matching on expression, breadth and length; a minimum of
100 matched-control draws; the matched-control band as the gate criterion; and
the stop-on-null rule that a result inside the band would not be expanded in
search of significance.

It does **not** designate macro-F1 as the gate metric, fix the treatment-set size
*k*, or set the final draw counts of 200 and 100. Those were decided after the
plan was written. `prespecification/results_addendum_2026-07-29.md` records every
such departure next to what was executed, and is the file to read second. The
plan itself is preserved unmodified, because a plan edited after seeing results
is not a plan.

The design was **not** lodged with an external registry, so the paper describes
it as *pre-specified* rather than pre-registered. A Git commit created now cannot
independently establish the 2026-07-07 date; the claim rests on the document
itself and on its use in the analysis scripts, not on repository history.

Note the contrast with the disease analysis: the specification hierarchy in
`E8_clinvar_adjusted.py` was fixed *after* exploratory results were seen, and
the script says so explicitly. Every specification is reported rather than the
primary being chosen blind.

---

## Reproducing

```bash
uv sync
```

Data acquisition (D-series) downloads model checkpoints and external databases
from their original sources. Model weights are not redistributed here.

Cost, on a single Apple Silicon workstation:

| Stage | Approximate time |
|---|---|
| Geometry screen, caller robustness, ESM-2 control | minutes to ~1 h |
| Covariate-aware annotation, diagnostics | seconds |
| **Matched deletion test (400 forward-pass sets)** | **~18 h** |

The deletion test is the only expensive step. It checkpoints every 10 control
sets and resumes on restart. `--primary-only` runs the treatment arm and its
matched null alone, skipping the sensitivity arm and k-sweep.

### Reproducibility notes

Two things bite in practice, and both are guarded in code:

- **Annotation table coverage.** `data/Table_S1.csv` at the repository root of
  the original project carried gene lengths for only 11,752 of 18,915 genes.
  Matching on median-imputed lengths made the controls *appear* balanced at
  SMD 0.070; the true figure for those controls was 0.334, and rebuilding them
  with complete lengths gave 0.077. Scripts now assert >=99% coverage and record the
  table's SHA-256 in every output.
- **Environment sensitivity.** The supervised probe baseline differs by 5.1e-4
  macro-F1 between CPU and MPS on the same machine and library stack. The pinned
  baseline is MPS, matching the device the ablations ran on, and reproduces
  bit-identically across two independent runs with the embedding cache cleared.
  BLAS threads are pinned and a full environment fingerprint is written
  alongside every result. The matched-null gate statistic compares treatment and
  control F1 directly and is invariant to this.

---

## Archived analyses

Superseded scripts and notebooks are retained under `notebooks/superseded/`
with a note explaining what replaced them. They are historical records, not
parts of the current analysis pipeline; `run_all.sh` does not execute them.

---

## Citation

```bibtex
@misc{whalley2026matchednull,
  author       = {Whalley, Justin P.},
  title        = {Comparative geometry of gene embeddings in single-cell
                  foundation models: model-specific outliers and limited
                  predictive value in Geneformer},
  year         = {2026},
  howpublished = {bioRxiv},
  doi          = {10.64898/2026.06.22.733850},
  note         = {Under consideration at Pacific Symposium on Biocomputing 2027}
}
```

Update to `@inproceedings` only on acceptance. See also `CITATION.cff`.

## License

See `LICENSE`. Model checkpoints and external databases retain their own
licences and are not redistributed here.
