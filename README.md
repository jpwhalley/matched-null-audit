# Matched-null auditing of gene-embedding outliers in single-cell foundation models

Code and reproducibility materials for:

> **A matched-null framework for auditing gene-embedding outliers in single-cell
> foundation models**
> Whalley, J.P. (2026). Submitted to *Pacific Symposium on Biocomputing 2027*.

Preprint: [doi:10.64898/2026.06.22.733850](https://doi.org/10.64898/2026.06.22.733850)
(v1 reported a different conclusion; see *Relationship to the earlier version*
below.)

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
> directly by `analysis/make_psb_figures.py` in the meantime. Outputs for the
> matched deletion test (E2, E7, E9 and Figure 4) are pending a corrected
> length-matched rerun and are deliberately absent rather than stale.

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

The pipeline runs left to right:

```
D-series  ->  analysis/E*.py  ->  outputs/*.{csv,json}  ->  notebooks/0*  ->  figures/
(acquire)     (compute)            (results)                (display)
```

---

## The five audit stages

The paper's contribution is the workflow, not the gene list. Each stage can
falsify a claim independently.

| Stage | Question | Implementation |
|---|---|---|
| 1. Geometry screen | Which genes are geometrically extreme? | `notebooks/P01_embedding_geometry.ipynb` |
| 2. Caller robustness | Is the outlier set an artefact of thresholding? | `analysis/E3_outlier_robustness.py` |
| 3. Non-transcriptomic control | Are the same genes outliers in protein-sequence space? | `analysis/E1_stage1_mapping.py`, `analysis/E1_stage2_esm2.py` |
| 4. Matched deletion test | Does removing them affect a downstream task more than matched controls? | `analysis/E2_downstream_ablation.py` |
| 5. Covariate-aware annotation | Does outlier status carry independent disease signal? | `analysis/E8_clinvar_adjusted.py`, `analysis/E6_class_association.py` |

Supporting diagnostics:

| Script | Purpose |
|---|---|
| `E7_cluster_metric_diagnostic.py` | Tests whether clustering metrics are usable as outcome measures. They are not; the paper reports macro-F1 alone. |
| `E9_token_occurrence_audit.py` | Measures actual token-level exposure of treatment vs control gene sets, to rule out an exposure-driven null. |

---

## Pre-specification

`prespecification/analysis_plan_2026-07-07.md` fixes the gate metric,
treatment-set size, number of matched-control draws, and the stop-on-null
decision rule for the deletion test (§2 line 33, §3 line 87). It predates all
deletion compute. It was **not** lodged with an external registry, so the paper
describes the design as *pre-specified* rather than pre-registered.

The file is preserved unmodified. Note that a Git commit created now cannot
independently establish its 2026-07-07 date; the claim rests on the document
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
- **Environment sensitivity.** The supervised probe baseline varies by
  ~6e-4 macro-F1 between execution environments while embeddings are unchanged.
  BLAS threads are pinned and a full environment fingerprint is written
  alongside every result. The matched-null gate statistic compares treatment and
  control F1 directly and is invariant to this.

---

## Relationship to the earlier version

The bioRxiv v1 preprint reported that embedding geometry predicts functional
fragility. That claim was tested here against expression-, breadth-, length- and
class-matched control gene sets and was not supported; it is retired. The
present work adds the non-transcriptomic control, the caller-robustness
analysis, the matched deletion test and the covariate-adjusted disease analysis.

Scripts and notebooks superseded by that revision are retained under
`notebooks/superseded/` with a note explaining what replaced them, rather than
deleted.

---

## Citation

```bibtex
@misc{whalley2026matchednull,
  author       = {Whalley, Justin P.},
  title        = {A matched-null framework for auditing gene-embedding outliers
                  in single-cell foundation models},
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
