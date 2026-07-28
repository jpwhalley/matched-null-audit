# Revision & Analysis Plan — "Glitch Genes" → Patterns

**Created:** 2026-07-07 (home computer, interactive session)
**Manuscript:** *Glitch genes: embedding geometry predicts functional fragility in single-cell foundation models* (Whalley, sole author)
**Trigger:** Genome Biology reject 2026-07-06 (2 reviewers; editor Zhicheng Ji). Transfer to BMC Artificial Intelligence offered — **declining**.
**Target venue:** **Patterns** (Cell Press). Backup: **TMLR**. Deliberately leaving the Springer Nature / NPG house after NM → NMI → GB.
**Deadline pressure:** **None.** Patterns has no submission window; the GB transfer clock is being declined; FPR AY26/27 credit already banked via the posted bioRxiv preprint (DOI 10.64898/2026.06.22.733850). Therefore: optimise for strength, not speed.

---

## 1. What the reviewers actually said (distilled)

Both reviewers called the idea creative and timely. Neither attacked analysis quality. They converged on **one inferential through-line**: the paper claims a geometric outlier is a *glitch* (representational failure) that *predicts functional fragility*, and the evidence supports "these are unusual, biologically central genes" but not "these are failures." Four load-bearing objections:

- **A. Extremeness ≠ failure (R1#2, R1#3).** No definition of what a "well-behaved" embedding space should look like, so an outlier can equally reflect genuine biological salience. Enrichment for constrained/disease/housekeeping/ribosomal/mito genes actually supports the *salience* reading — cutting against the failure framing.
- **B. No downstream link (R1#4, R2#3).** Title says "predicts functional fragility" but there's no demonstration that geometry degrades a real task (annotation, clustering, perturbation) beyond an expression-matched control, or that it can't be fixed by calibration.
- **C. The paper's own control undercuts the claim (R2 minor#2).** The exposure regression shows the anomaly–perturbation correlation *largely disappears* after controlling for expression level, breadth, and gene length. Reviewer read: geometry may be a proxy for training-exposure imbalance, not an independent failure mode.
- **D. The glitch analogy is inverted and imported uncritically (R1#1, R2#3/4).** NLP glitch tokens come from *under*-training; these outliers are *high*-exposure genes. Same word, opposite mechanism → the borrowed hubness/anisotropy→fragility intuition doesn't transfer for free.

Constructive/addressable methodology points:
- **E. z>3 is fragile on non-Gaussian distributions (R2#5)** — Geneformer norm right-skewed, scGPT norm bimodal. Wants MAD / percentile / mixture-model / empirical-null alternatives.
- **F. Panel is dated (R2#1, R2#4)** — wants newer / architecturally-diverse / LM-initialised models.
- **G. ESM-2 sequence-embedding control (R2#2)** — the sharpest suggestion; see E1.
- **H. Synthetic-cell methods underspecified (R2 minor#1).**

## 2. The claim-level decision (the fork everything hangs on)

The revision must pick which paper this is, and the experiments in §3 decide it:

- **"Predictor" paper** — keeps a causal/predictive claim. *Only tenable if* the E2 pilot shows outlier genes degrade a real task beyond controls matched on expression, breadth, length **and gene class**. Then "predicts functional fragility" is earned.
- **"Diagnostic / characterization" paper — treat as the likely landing, and an equally strong paper.** Geometric outliers in scFM embeddings largely track training-exposure imbalance and intrinsic biological salience; auditing geometry is a diagnostic lens that repeatedly surfaces biologically special genes, not a fragility predictor. This is honest, publishable, and a *better* Patterns/TMLR fit than an over-claimed predictor paper. **We are not "rescuing" a predictor paper; we are letting the evidence choose between two good papers.**

**Pre-commitment (non-negotiable):** if the E2 pilot effect is not outside the matched-control null, **stop defending "predicts fragility" and pivot cleanly** — do not escalate the pilot into a giant undifferentiated experiment fishing for significance. R2 already found the thread (C); a buried null becomes reject #4.

**On the inversion (D) — turn the criticism into the thesis.** The paper currently both *leans on* the glitch-token analogy to import the fragility intuition *and* celebrates that glitch genes are its mechanistic opposite. Can't do both. Resolution: **foreground the inversion as a positive finding** ("we expected NLP-style under-training glitches; we found the opposite — high-exposure genes are the geometric outliers"), and make the fragility claim (if retained) stand entirely on E2's independent evidence, not on the analogy. This directly answers D and repurposes R1#1.

**Terminology / retitle — DECIDED (2026-07-07): retire "Glitch Genes" as the paper's banner now.** The reviewers are telling us the term is doing active damage (imports the inverted intuition it can't support). Keep it only as (a) private historical shorthand and (b) at most one sentence in the intro noting the NLP lineage that motivated the audit. Main title from §6.

## 3. The experiments (specified to run)

Inputs already on disk (READ-ONLY `coding/a_bio_AI/data/`): `gene_embeddings.npy` (Geneformer), `scgpt_gene_embeddings.npy`, scFoundation equivalents (from `sf_*` notebooks), `gene_names.json`, `outlier_genes.json`, `gene_embedding_geometry.csv` (+ scgpt/sf variants), `exposure_regression_results.csv`, `perturbation_sensitivity.csv`, `pbmc3k_raw.h5ad`. **New outputs must go to a writable location** (see §7 blocker).

### Execution order (staged, gate-driven — cheapest & most protective first)

The experiments are **not** run in parallel or in E1→E5 label order. Each stage can change or stop the next:

1. **E3 robustness FIRST.** Cheapest, and it protects everything downstream. If the outlier set is unstable under MAD / percentile / mixture thresholds, the whole revision must soften *before* we spend compute on E1/E2. If stable, it's a strong standing reply to R2#5.
2. **E1 Stage 1 — mapping feasibility** (metadata only, no embeddings). Decide whether the ESM-2 comparison is even clean/unbiased enough to run.
3. **E1 Stage 2 — ESM-2 embeddings + audit**, only if Stage 1 is clean enough.
4. **E2 pilot — the claim gate.** Tight pilot (1–2 models, PBMC3k + one external labelled set). Result decides predictor-vs-characterization framing. **Stop-on-null** rule applies.
5. **Decide title/framing from E1 + E2 jointly.**
6. **E5 (newer model) — held.** Only if steps 1–4 are strong *and* the paper still needs breadth. The reviewers' central objection is "you haven't shown these outliers matter," not "only three models" — fix that first; don't let model-zoo expansion become a swamp.

**Gate discipline — one dated verdict memo per gate.** After E3, after E1-mapping, and after the E2 pilot, write a ≤1-page dated memo (`E3_verdict_<date>.md`, `E1_mapping_verdict_<date>.md`, `E2_pilot_verdict_<date>.md`) in the writable folder, each ending with an explicit **CONTINUE / PIVOT / STOP** verdict and one line of rationale. This is the structural guard against the revision turning into an emotional salvage operation: the decision is written down and dated before the next stage's compute is spent.

### E1 — ESM-2 sequence-embedding control *(answers G, strengthens A/C/D) — TWO STAGES*
**Question:** are the same genes flagged as outliers by an embedding space that never saw expression data?

**Stage 1 — mapping feasibility / metadata pass (NO embeddings yet).** Cheap, and it decides whether Stage 2 is even worth running unbiased. See **Appendix A** for the concrete recipe. Deliverable = a feasibility table answering:
- What fraction of each model's vocab maps cleanly to a canonical protein? Overall vs **within the outlier set vs non-outlier set** (the key bias check).
- Which categories drop out? Non-coding genes have *no* protein and can never appear in ESM-2 space — if outliers are disproportionately non-coding, that's a structural scope limit to state, not hide.
- Is the comparison biased toward protein-coding / high-confidence genes? If so, restrict E1 to the shared mappable protein-coding subset and declare that scope explicitly.
- **Gate:** proceed to Stage 2 only if mapping is clean enough (target: majority of protein-coding outliers map, coding/non-coding split not wildly skewed between outlier/non-outlier sets).

**Stage 2 — ESM-2 embeddings + audit (only if Stage 1 clean).**
1. Generate per-protein ESM-2 embeddings (esm2_t33_650M_UR50D; mean-pool residues, exclude BOS/EOS). Cache once. Handle >1022-residue proteins by windowed mean-pool (Appendix A).
2. Run the *identical* four-metric audit + same outlier caller on the ESM-2 space, restricted to the shared mappable subset.
3. Compare vs each scFM: outlier-set overlap (Jaccard + hypergeometric), Spearman rank-correlation of composite anomaly scores, and category concordance (do housekeeping/ribo/mito/constrained surface again?).
**Interpretation (both outcomes are findings, neither kills the paper):**
- *Divergent outliers* → outliers are model/tokenization-specific → supports the audit/fragility story. **Best case.**
- *Same outliers* → they reflect intrinsic biological properties (length, composition, conservation) repeatedly pushed to representational extremes → **decisively reframes the thesis to "geometry reveals biologically special genes," not "model failure."** This is a strong, honest result — arguably the more interesting one.
**Effort:** Stage 1 low (metadata only); Stage 2 moderate (ESM-2 inference). **Predict the direction with Justin before Stage 2** (§8).

### E2 — Expression-matched downstream ablation *(THE claim gate; answers B and C) — PILOT FIRST*
**Question:** does removing/masking high-anomaly genes degrade a real task *beyond* controls matched on the confounds that killed the correlation in C — **and** beyond gene-class confounds?

**Run as a tight pilot, not a giant experiment.** Scope the pilot deliberately small so it can cleanly gate the framing:
1. **Models:** 1–2 (start with the model whose outlier set is most stable per E3; likely Geneformer + one).
2. **Data:** PBMC3k **plus one external labelled dataset** (e.g. a Tabula Sapiens subset or an annotated immune atlas) — enough to show it's not PBMC-specific, no more.
3. **Treatment set** = top-k high-anomaly genes. **Ablate** identically (mask / delete-and-reindex per tokenization) for treatment and each control draw.
4. **Matched controls — three-covariate baseline + a gene-class layer (this is the key hardening).** Reviewers will note that many outliers are housekeeping / ribosomal / mitochondrial / constrained genes, so matching on expression, breadth, and length alone is dismissable. Add one of:
   - **Primary: class-stratified matching** — match each outlier to controls *within the same gene class* (housekeeping, ribosomal, mito, constrained, other), so the comparison is outlier-vs-non-outlier *within* class.
   - **Sensitivity: drop ribosomal + mitochondrial outliers and re-test** — if the effect survives with the "obvious" categories removed, it's much harder to dismiss.
   Use ≥100 bootstrap matched draws for a null band.
5. **Tasks (≥2):** (a) cell-type annotation accuracy (zero-shot / linear-probe on cell embeddings); (b) clustering stability (ARI / NMI vs reference labels). Report effect size + CI, not just p.
**Gate criterion:** treatment Δ **outside** the matched-control null band (and surviving the class layer) → independent fragility signal survives both confounds → "predictor" framing earned → *then* scale to full model set.
**Stop-on-null (pre-committed):** treatment Δ **inside** the band, or effect evaporating once ribo/mito are removed → geometry adds nothing beyond expression + class → **pivot cleanly to characterization framing (§2); do not expand the pilot chasing significance.**
**Effort:** moderate as a pilot (the whole point of piloting first); only scales to "high" if the gate opens.

### E3 — Non-parametric outlier robustness *(RUN FIRST; cheap, protective; answers E)*
Re-derive outliers with MAD-based z, percentile thresholds, Gaussian-mixture / empirical-null calibration per metric. Report outlier-set stability (Jaccard vs the |z|>3 set) and re-run the headline enrichments on each. Show conclusions don't hinge on the Gaussian assumption on skewed/bimodal metrics. **Gate:** if the outlier set is unstable across methods, the entire revision must soften its claims *before* E1/E2 compute is spent; if stable, this is a standing reply to R2#5 and the anchor set for E1/E2. **Low effort, high payoff — do it first.**

### E4 — Synthetic-cell methods + real-data consistency *(answers H)*
Document the synthetic-cell generation exactly (rank-value construction, expression/rank distribution preservation for Geneformer's rank encoding), and add a side-by-side that synthetic-cell perturbation results track real single-cell results (already have `perturbation_synthetic_vs_real.csv`). Mostly writing + one confirmatory figure.

### E5 — Additional / newer model *(HELD — do not run unless E1–E3 strong AND breadth still needed; answers F)*
Three architectures already demonstrate cross-architecture generality, and an AI/methods audience discounts "not the newest model" more than biologists do. The reviewers' central objection is "you haven't shown these outliers matter," **not** "only three models" — so fix E2 first. New models can become a swamp. **Do not run pre-emptively.** If steps 1–4 land strong and the paper genuinely needs breadth, pre-scope *one* drop-in candidate (a 2025 scFM or an LM-initialised model, which also answers R2#4's "LM-derived model" ask) — not a zoo.

## 4. Narrative / framing changes (independent of experiment outcomes)

- **De-glitch** the title and text; keep a one-line nod to the NLP analogy only to *set up the inversion as a finding*.
- **Define the null model** ("what a well-behaved embedding space should look like") explicitly up front — answers A directly. E1 and E3 supply the empirical baselines.
- **Promote the expression confound (C) from a buried control to a central, honest result.** Frame it as: geometry is substantially — but test whether *entirely* — explained by exposure; E2 is the arbiter.
- **Reframe the developer recommendation.** Drop "include geometric outlier reports as standard supplementary material" unless E2 links geometry to a task metric; otherwise present the audit as an *interpretability/diagnostic lens*, not a QC gate.

## 5. Reviewer-response coverage map

| Reviewer point | Addressed by |
|---|---|
| R1#1 inverted analogy | §2 reframe (inversion-as-finding), §4 de-glitch |
| R1#2 no "well-behaved" baseline | §4 null model + E1 + E3 |
| R1#3 salience vs failure | E1 (intrinsic-vs-model test) + E2 |
| R1#4 no downstream/calibration link | **E2** |
| R2#1 dated panel | §3 E5 (optional) + cross-architecture argument |
| R2#2 ESM-2 control | **E1** |
| R2#3 downstream examples | **E2** |
| R2#4 LM-derived model | E5 candidate = LM-initialised model |
| R2#5 z>3 fragility | **E3** |
| R2 minor#1 synthetic cells | E4 |
| R2 minor#2 expression confound | §4 (promote to central) + E2 arbitration |

## 6. Title candidates (drop "Glitch Genes")

- *Auditing the geometry of gene embeddings in single-cell foundation models*
- *Geometric outliers in single-cell foundation model embeddings reflect training exposure, not representational failure* (if E2 null)
- *When does embedding geometry predict fragility? A weight-only audit of single-cell foundation models* (if E2 positive)
- *A weight-only geometric audit of gene representations across single-cell foundation models*

## 7. Timing & sequencing (recommended)

**Reject the rush.** No deadline; banked credit; strength wins.

1. **This week (home):** this plan + Appendix A recipe + reframe/title drafting. *(No compute, no write access needed.)*
2. **From Jul 17 (office — write access + real compute), in staged order:** E3 (robustness, first) → E1 Stage 1 (mapping feasibility) → E1 Stage 2 (ESM-2, if clean) → E2 pilot (claim gate, stop-on-null) → E4 (synthetic-cell writeup, alongside). Each stage's result can stop or reshape the next; do not batch them.
3. **Late Jul → mid-Aug (perturb_ai MLCB reviews land; decision Aug 14):** use E1+E2 results **plus** that independent second reviewer set to lock claim level (predictor vs characterization), finalise framing/title, write cover letter, submit to Patterns. E5 only if still needed.

**Blocker requiring Justin's decision:** `coding/a_bio_AI` is READ-ONLY per access rules. To run E1–E4 and save outputs, set up a **writable exploration folder** (mirror pattern, e.g. `coding/a_bio_AI_revision/` or under `writing/manuscripts/bio_ai/analysis/`) that can read the existing embeddings/data. Nothing gets written into the read-only tree.

## 8. Open decisions for Justin

1. **ESM-2 expected direction?** Same outliers (→ characterization) or different (→ model-specific fragility)? Worth predicting before Stage 2. *(Still open — Justin to call.)*
2. **Writable folder** — approve location for experiment outputs (§7 blocker). *(Still open — required before Jul 17 execution.)*
3. ~~**Retitle**~~ — **RESOLVED 2026-07-07:** retire "Glitch Genes" as banner now (§2, §6).
4. ~~**E5**~~ — **RESOLVED 2026-07-07:** held; do not run pre-emptively (§3 E5).
5. **Merge-with-perturb_ai contingency** — confirmed last-resort drawer item, not a plan (lineage doc warns of pattern-match discounting).

---

## Appendix A — ESM-2 gene→protein mapping recipe (E1 Stage 1, then Stage 2)

Goal: a reproducible, versioned map from each model's gene vocabulary to a single canonical protein sequence, with an explicit bias audit — *before* spending any ESM-2 inference.

### A.0 Inputs (read-only)
- Per-model vocab + outlier calls: `coding/a_bio_AI/data/gene_names.json`, `outlier_genes.json`, `gene_embedding_geometry.csv` (+ `scgpt_*`, `sf_*` variants).
- Identifier types differ by model — **normalise per model, don't assume a shared ID space:**
  - **Geneformer** — Ensembl gene IDs (ENSG). Strip version suffixes (`ENSG00000123456.7` → `.7` dropped) before mapping.
  - **scGPT** — HGNC gene symbols; resolve aliases/withdrawn symbols to current HGNC before mapping.
  - **scFoundation** — gene symbols (~19,264-gene panel); same alias normalisation.

### A.1 Canonical protein assignment (one protein per gene)
Use **MANE Select v1.4** as the primary canonical source (NCBI/EMBL-EBI agreed RefSeq+Ensembl transcript per gene; versioned → reproducible). `MANE.GRCh38.v1.4.summary.txt.gz` gives Ensembl gene ↔ Ensembl protein (ENSP) ↔ RefSeq protein (NP_) ↔ symbol.
- Geneformer: ENSG → MANE row → canonical protein.
- scGPT / scFoundation: symbol → HGNC → ENSG → MANE row.
- **Fallback** for genes absent from MANE (rare, or non-MANE canonical): UniProtKB **reviewed (Swiss-Prot)** canonical isoform via UniProt ID-mapping (Ensembl_gene/Gene_Name → UniProtKB-Swiss-Prot); take the single reviewed canonical entry.
- Sequence source: the MANE protein sequence (or the UniProt canonical FASTA for fallbacks). Keep source consistent per gene and record which was used.
- Tooling: `mygene` (mygene.info) or `pybiomart` for ID resolution over the sandbox's allowlisted network; or a **static cached MANE summary + UniProt idmapping table** if network is restricted (preferred for reproducibility — download once, commit to the writable folder). Do NOT use raw curl/wget per harness rules; use the packages' HTTP clients or a pre-cached file.

### A.2 Feasibility table (the Stage-1 deliverable)
One row per (model, gene) with columns:
`model, vocab_id, gene_symbol, ensembl_gene, biotype, mane_status, mapped_protein_id, protein_len, in_outlier_set, gene_class, mapping_status`
where `gene_class ∈ {housekeeping, ribosomal, mitochondrial, constrained, transcription_factor, other}` (housekeeping = HRT-Atlas / Eisenberg–Levanon; ribosomal = RPL/RPS + MRPL/MRPS; mito = MT- + MitoCarta nuclear-encoded; constrained = gnomAD pLI>0.9 or LOEUF<0.35).

### A.3 Bias audit (decides the Stage-2 gate)
- Mapping rate: overall, **and split outlier vs non-outlier**, **and per gene_class**.
- Coding vs non-coding × outlier vs non-outlier contingency (Fisher/χ², report OR) — the structural bias check. Non-coding outliers cannot enter ESM-2 space; quantify how much of the outlier set that removes.
- **Gate:** proceed to Stage 2 only if the majority of *protein-coding* outliers map and the coding/non-coding split isn't wildly skewed between sets. Whatever the shared mappable subset is, **that** becomes the declared scope of E1; state it plainly.

### A.4 Stage 2 — embeddings (only past the gate)
- Model: `facebook/esm2_t33_650M_UR50D` (650M, 1280-d) via HuggingFace or fair-esm. (Note `esm2_t36_3B` as an optional robustness check, not default.)
- Per protein: forward pass, **mean-pool residue representations excluding BOS/EOS** → one 1280-d vector. Cache to `.npy` keyed by vocab_id.
- Long proteins (>1022 residues, ESM-2 context limit): **windowed mean-pool** (overlapping 1022-residue windows, mean each, then average) rather than truncation, to avoid biasing long proteins (many outliers are long/structural).
- Then run the identical four-metric audit + outlier caller on this matrix, restricted to the shared mappable subset, and compute the A.3 comparisons (Jaccard, hypergeometric, Spearman of composite scores, category concordance) against each scFM's outlier set.
