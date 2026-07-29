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

Both limbs were met. The pre-committed pivot was executed.

| | full (50 genes) | sensitivity (36 genes) |
|---|---:|---:|
| treatment − control mean | **−0.002946** | **+0.000544** |
| gate statistic z | **−1.5619** | **+0.3481** |
| controls at least as damaging | **11/200** | **132/200** |
| one-sided empirical p (add-one) | **0.0597** | **0.6617** |
| inside 95% band | yes | yes |

Verdict: **NULL**. The treatment did not exceed the matched-control null on
either arm.

---

## Fixed in advance

- class-stratified matching on expression, breadth and length
- ≥100 matched draws (200 were run)
- gate criterion: treatment Δ outside the matched-control band
- the two-limbed stop-on-null rule above

## Decided after the plan was written

Recorded here rather than presented as pre-specified.

| Choice | Plan | Executed |
|---|---|---|
| Outcome metric | "Tasks (≥2): (a) annotation accuracy; (b) clustering stability (ARI/NMI)" — no gate designated | retrained macro-F1 designated as the gate |
| Treatment size *k* | "top-k", k free | k = 50, with k = 25 and 100 as a descriptive sweep |
| Draw count | ≥100 | 200 |
| External dataset | "PBMC3k **plus one external labelled set**" | PBMC3k only; the Tabula Sapiens attempt is described in the manuscript's limitations |
| MANE release | v1.4 (Appendix A.1) | v1.3 |

Clustering was a **pre-specified outcome** that is excluded from the headline
figure. The exclusion rests on precision — the ARI null band spans ~42% of its
baseline against ~0.8% for macro-F1 — and is reported rather than applied
silently. An earlier justification, that random matched deletion reliably
*improves* ARI, did not survive the corrected controls and is not relied on.

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

---

## Not done, deliberately

The empirical p of 0.0597 sits near 0.05. **No additional control draws were
generated after seeing it.** The pre-specified criterion is the 2.5th percentile
of the matched-control distribution and it was met; extending the null after
observing a borderline p would be outcome-driven whichever way it moved.
