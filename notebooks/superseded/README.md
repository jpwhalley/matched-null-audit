# Superseded notebooks

These reproduce analyses from the bioRxiv v1 preprint that the revision either
retired or replaced. They are kept rather than deleted so that a reader can see
what changed and why, and can re-run the earlier analysis if they wish to check
the correction for themselves.

**None of these produce results reported in the current manuscript.**

| Notebook | Why it was superseded | What replaced it |
|---|---|---|
| `05_perturbation_sensitivity.ipynb`<br>`P06_perturbation_sensitivity.ipynb` | Reported ρ = 0.725 between geometric anomaly and in-silico perturbation sensitivity, and treated it as evidence that geometry predicts functional fragility. The correlation is real but was never tested against controls matched on expression, breadth, length and gene class. | `analysis/E2_downstream_ablation.py` — the matched deletion test. The effect falls inside the matched-control null. |
| `07_exposure_regression.ipynb` | Showed the anomaly–perturbation partial correlation collapsing to −0.019 once expression was controlled, and reported ΔR² = 0.004 from adding anomaly score to an expression-only model. Presented as a limitation. | Not replaced — **promoted**. This result is now the organising argument of the paper, and its logic is carried by the matched-null design throughout. |
| `02_metric_decomposition.ipynb` | Separated magnitude-driven from isolation-driven outliers and characterised each. Sound, but framed as identifying two "failure modes" when no failure had been demonstrated. | Retained as supplementary characterisation; the "failure mode" framing is dropped. |
| `04_biological_annotation.ipynb` | Reported unadjusted class enrichments, including disease association, without covariate adjustment or a stated class scheme. | `analysis/E8_clinvar_adjusted.py` (covariate-adjusted, Firth primary) and `analysis/E6_class_association.py` (both class schemes reported). |
| `08_external_validation.ipynb` | Combined gnomAD, DepMap, IMPC, ClinVar and Replogle results as "external validation". Mixed significant enrichment, null results and a correlation with a different biological quantity under one heading, all unadjusted. | Not carried forward. A covariate-aware reanalysis would be required before any of it could support a claim. |

## The short version

The v1 analysis established that geometric outliers exist and are biologically
non-random. It then inferred that they are functionally fragile. That inference
was the weak step, and it is the one the revision tested and rejected.

The exposure regression in `07` had already shown the confound clearly. In v1 it
appeared as a caveat. It is now the point.
