# Results addendum — 2026-07-29

**The dated plan `analysis_plan_2026-07-07.md` is preserved unmodified.** This
addendum records what was executed against it and where execution departed from
it. Nothing here is a revision of the plan; a plan edited after seeing results
is not a plan.

---

## Outcome under the pre-specified rule

The stop-on-null rule has two limbs, joined by *or*:

> treatment Δ **inside** the band, **or** effect evaporating once ribo/mito are
> removed → geometry adds nothing beyond expression + class → pivot cleanly to
> characterisation framing; do not expand the pilot chasing significance.

Both limbs were met and the pre-specified pivot was executed.

| | full (50 genes) | sensitivity (36 genes) |
|---|---:|---:|
| treatment − control mean | **−0.002946** | **+0.000544** |
| gate statistic z | **−1.5619** | **+0.3481** |
| controls at least as damaging | **11/200** | **132/200** |
| one-sided empirical p (add-one) | **0.0597** | **0.6617** |
| inside 95% band | yes | yes |

Verdict: **NULL**. The treatment did not exceed the matched-control null on
either arm.

### Replication (Tabula Sapiens Immune, primary arm only)

| | full (50 genes) |
|---|---:|
| baseline macro-F1 | 0.6348 |
| treatment − control mean | **−0.007894** |
| control SD | 0.006284 |
| gate statistic z | **−1.2561** |
| controls at least as damaging | **11/100** |
| one-sided empirical p (add-one) | **0.1188** |
| inside 95% band | yes |

Concordant gate decision on a lower-baseline task. The raw contrast is larger
than PBMC3k's, inside a null that is intrinsically more variable, so it is less
exceptional after standardisation. Executed on separate hardware under Python
3.11.14. The sensitivity arm was not replicated.

---

## Fixed in advance

- class-stratified matching on expression, breadth and length
- ≥100 matched draws (PBMC3k: 200; Tabula Sapiens: 100)
- gate criterion: treatment Δ outside the matched-control band
- the two-limbed stop-on-null rule above

## Decided after the plan was written

Recorded here rather than presented as pre-specified.

| Choice | Plan | Executed |
|---|---|---|
| Outcome metric | "Tasks (≥2): (a) annotation accuracy; (b) clustering stability (ARI/NMI)" — no gate designated | retrained macro-F1 designated as the gate |
| Treatment size *k* | "top-k", k free | k = 50, with k = 25 and 100 as a descriptive sweep |
| Draw count | ≥100 | PBMC3k 200; Tabula Sapiens 100 |
| External dataset | "PBMC3k **plus one external labelled set**" | Met in part: PBMC3k plus a Tabula Sapiens Immune subset (22 of 43 types, 3,919 cells), primary arm only |
| MANE release | v1.4 (Appendix A.1) | v1.3 |

Clustering was a **pre-specified outcome** that is excluded from the headline
figure. The exclusion rests on precision — the ARI null band spans ~42% of its
baseline against ~0.8% for macro-F1 — and is reported rather than applied
silently. An earlier justification, that random matched deletion reliably
*improves* ARI, did not survive the corrected controls and is not relied on.

On Tabula Sapiens, cluster NMI *does* meet the precision criterion (band 3.9% of
baseline) where it failed on PBMC3k (10.6%), but shows no treatment effect
against its controls (z = −0.06). It remains secondary because macro-F1 was
fixed as the gate before Tabula Sapiens was analysed. Reported rather than acted
on: changing the outcome metric after seeing a second dataset is the move this
addendum exists to make visible.

---

## Corrections applied after the plan

1. **Gene-length coverage.** The annotation table first used carried lengths for
   11,752 of 18,915 genes, median-imputed downstream, which silently broke the
   length-matched design. Controls that *appeared* balanced at SMD 0.070 scored
   0.334 against complete lengths; rebuilt with complete lengths they return
   0.077. Scripts now assert ≥99% coverage and record the table's SHA-256.

   **This correction moved the primary result toward the boundary, not away from
   it.** Better-matched controls do less damage on average, so the treatment
   stands out more: z went from −1.41 to −1.562 and the empirical p from 0.080 to
   0.0597. The sensitivity arm's margin roughly halved. Recorded because the
   honest form of "we corrected our controls" includes what the correction cost.

2. **Baseline provenance.** The probe baseline differs by 5.1e-04 macro-F1
   between CPU and MPS on the same machine and library stack. The pinned
   baseline is MPS, matching the device the ablations ran on, and reproduces
   bit-identically across two independent runs with the embedding cache cleared.
   An earlier baseline value of 0.9242120 carried no environment fingerprint and
   could not be reproduced on either device; its provenance is unrecoverable and
   it is not used. The gate statistic, rank and empirical p are invariant to a
   common baseline shift and are unaffected.

3. **Empirical p definition.** Reported as (b+1)/(m+1), the add-one form, since
   the controls are a random sample of the possible matched sets.

4. **Ribosomal panel definition.** The panel was defined by the pattern
   `^(RPL|RPS|MRPL|MRPS)\d`, which requires a digit immediately after the
   family prefix and therefore misses four genuine ribosomal protein genes:
   RPLP0, RPLP1, RPLP2 and RPSA. All four are Geneformer geometric outliers.
   The pattern was corrected on 2026-07-30 in E1, E3, E6 and E8.

   **This does affect the deletion controls, and the first version of this
   note was wrong to say otherwise.** No treatment gene changes class, so the
   treatment strata are unchanged. But three of the four (RPLP0, RPLP2, RPSA)
   were eligible as `constrained` controls and one (RPLP1) as `other`, and they
   were drawn repeatedly: 314 occurrences across the 200 PBMC3k full-treatment
   control sets, 325 across the sensitivity sets, and 29 across the 100 Tabula
   Sapiens sets. Under the corrected pattern they compete only for the 13
   ribosomal slots, so the eligible pool and therefore the null distributions
   change. The net effect on balance is not
   predictable in advance: the four leave the constrained and other pools but
   enter the ribosomal pool, so matching may improve, worsen or be unaffected on
   any given covariate. It will be measured, not assumed.

   The reported deletion results were produced under the narrower pattern.
   Until the controls are rebuilt and the ablations rerun, the manuscript states
   this explicitly in the matched-control Methods rather than implying a single
   definition throughout.

   Effect on the analyses that were rerun (E1, E3, E6, E8), all small and none
   changing a verdict. The deletion test is not included here and its exposure
   to the correction is unresolved:
   ribosomal enrichment OR 36.7 to 39.3; constrained (mutually exclusive) 1.59
   to 1.53; adjusted ClinVar OR 1.05 to 1.06 (Geneformer), 1.22 to 1.23
   (scGPT), 0.98 unchanged (scFoundation); exploratory shared set 1.45 to 1.47.
   E3 remains MIXED and E1 remains DIVERGENT.

---

## Pre-commitments for the ribosomal rerun, recorded 2026-07-30 before any control was drawn

The correction in item 4 changes the eligible control pool, so the matched nulls
are being rebuilt. Because the current result is already known, the terms are
fixed here and committed before the rerun starts. This section is written in
advance and is not to be edited afterwards.

1. **The pinned HGNC panel is primary.** Ribosomal membership is defined by
   `data/ribosomal_panel.csv`, built from HGNC curated gene groups 728, 729 and
   646, retrieved 2026-07-30, SHA-256 `acd2a7ad...5cf5842e`, 171 genes, all 171
   present in `Table_S1` and agreeing on Ensembl and symbol keys. Membership is
   resolved by Ensembl gene ID. This replaces the symbol regex entirely, in E2
   as well as E1/E3/E6/E8. Relative to the interim four-gene regex fix it adds
   AURKAIP1, CHCHD1, DAP3, FAU, GADD45GIP1, PTCD3 and UBA52, and removes ten
   RPS6K* kinases, RPS19BP1 and the obsolete symbol MRPS36. The treatment set
   is unchanged: no panel change falls inside the top 50, and UBA52 and FAU are
   outliers ranking 119 and 304.
2. **Nothing else changes.** The gate metric (retrained macro-F1), the treatment
   sets, the matcher, the random seed and the draw counts (PBMC3k 200 full and
   200 sensitivity, Tabula Sapiens 100 primary) are unchanged. Baselines and
   treatment selection are not recomputed.
3. **The corrected results replace the current ones regardless of direction or
   significance.** If the treatment moves outside the matched-control band, the
   paper reports that, with the same restraint it currently applies to the null.
4. **The digit-anchored results are archived, not discarded.** They remain in
   the repository as a disclosed provenance and sensitivity analysis.
5. **No additional draws will be generated after the corrected result is seen.**
   The criterion remains the 2.5th percentile of the matched-control
   distribution.

`E2_downstream_ablation.py` records the panel SHA-256 beside every control
cache and refuses a cache drawn under a different panel or draw count; ablation
checkpoints are bound to the SHA-256 of the control set they were written
against. The old and new nulls cannot be mixed or spliced.

Why a curated panel rather than a better regex: a regex on symbols does not
track curation. The digit-anchored form missed four genuine members; correcting
it to admit them still admitted ten kinases and an obsolete symbol. The failure
was the method, not the expression, so the expression was replaced.

---

## Matching balance under the HGNC panel, recorded 2026-07-30 before any ablation completed

Rebuilding the PBMC3k nulls under the pinned panel moved expression balance in
the wrong direction, by a modest amount:

| arm | v1 (digit-anchored) | v2 (HGNC panel) |
|---|---:|---:|
| PBMC3k full, 50 genes | SMD 0.2881 | SMD 0.3152 |
| PBMC3k sensitivity, 36 genes | SMD 0.3585 | SMD 0.3982 |
| Tabula Sapiens full, 50 genes | SMD 0.2828 | SMD 0.283 |
| Tabula Sapiens sensitivity, 36 genes | SMD 0.3097 | SMD 0.313 |

Class strata remain exact throughout, and breadth and length are unaffected
(PBMC3k full: SMD 0.021 and 0.077; Tabula Sapiens full: 0.024 and 0.112).

The degradation is specific to PBMC3k. Tabula Sapiens moved by less than 0.004
on both arms, because the two datasets draw from different candidate pools
after the gene-length filter (18,333 eligible for PBMC3k against 18,909 for
Tabula Sapiens). Whatever the panel change costs in matching quality, it is not
a general property of the panel.

The mechanism is not the obvious one. Removing the ten RPS6K kinases cleans the
ribosomal stratum, but those twelve removals move into the constrained, disease
and other pools while the seven additions leave them. Those pools supply 19
constrained and 9 other control slots per draw, so they lost their most highly
expressed members and gained moderately expressed ones. Mean control expression
fell from 3.92 to 3.74 on the full arm and from 2.04 to 1.79 on the sensitivity
arm.

**Direction of the resulting bias, stated before the result is known.** Controls
now sit further below the treatment on expression. Lower-expression genes
contribute fewer tokens per cell, so matched controls should damage annotation
less than before, which makes the treatment appear relatively more damaging and
drives the gate statistic more negative. The imbalance therefore biases *toward*
exceeding the matched-control null, not toward it.

Consequently: a null result under the corrected panel is stronger evidence than
the v1 null, because it survives worse expression matching. A result that
exceeds the null must be reported together with this imbalance as a candidate
explanation, and not presented as clean evidence of a geometry effect. This
paragraph is written while the ablation is still running and is not to be
revised after the outcome is seen.

---

## Not done, deliberately

The empirical p of 0.0597 sits near 0.05. **No additional control draws were
generated after seeing it.** The pre-specified criterion is the 2.5th percentile
of the matched-control distribution and it was met; extending the null after
observing a borderline p would be outcome-driven whichever way it moved.

---

## Outcome under the corrected HGNC panel (recorded 2026-07-31)

Both datasets were re-run against matched nulls rebuilt from the pinned HGNC
panel (groups 728, 729, 646; 171 genes; sha256 `acd2a7ad…`). The pre-specified
rule was applied unchanged.

| arm | treat − ctrl mean | z | at least as damaging | add-one p | band | verdict |
|---|---|---|---|---|---|---|
| PBMC3k primary | −0.002805 | −1.500 | 12 / 200 | 0.0647 | [−0.003586, +0.003022] | inside |
| PBMC3k sensitivity | +0.001216 | +0.574 | 145 / 200 | 0.7264 | [−0.007645, +0.001089] | inside |
| Tabula Sapiens primary | −0.007297 | −1.120 | 14 / 100 | 0.1485 | [−0.012209, +0.011704] | inside |

Gate: NULL on every arm. The Tabula Sapiens sensitivity controls were drawn
(100 draws of 36 genes) but the arm was not evaluated; the replication is
primary-only and is labelled as such in the manuscript and in Figure 4.

The prospectively recorded directional prediction did not hold for the
standardised contrast: z moved from −1.562 to −1.500 and the add-one empirical
p-value from 0.0597 to 0.0647. However, the absolute margin above the
pre-specified 2.5th-percentile boundary narrowed from 0.00104 to 0.00044. The
treatment remained inside the matched null under both specifications.

The narrowing is driven by the null rather than the treatment: the treatment
moved 0.00005, while the 2.5th percentile of the control distribution rose
from −0.00419 to −0.00359. The rebuilt controls contain fewer extremely
damaging draws. This does not establish that the proposed mechanism is wrong;
it shows that token count alone did not predict the net macro-F1 effect once
control composition changed.

Under pre-commitment 5, no additional draws were generated after these results
were seen. The v1 digit-anchored controls are retained under
`archive/ribo_v1_digit_anchored/` and are labelled superseded, not current.

**Scope of this record.** This comparison is provenance. The manuscript
reports the corrected result and its bounded conclusion only; the bioRxiv
changelog states that the corrected control definition changed the numerical
estimates but not the verdict. Neither carries the v1 statistics, the boundary
margin or the directional prediction.
