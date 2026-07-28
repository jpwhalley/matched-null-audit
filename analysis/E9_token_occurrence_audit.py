"""
E9 — Token-occurrence audit for the matched-deletion design.

WHY THIS EXISTS
---------------
Matching on global expression level, breadth and length does not guarantee that
treatment and control genes are equally PRESENT in the cells actually analysed.
A gene only influences a prediction if its token appears in that cell's
rank-value encoding. If the treatment set were systematically less represented
in the input than its controls, the deletion test would remove less signal and
the observed null could be an exposure artefact rather than a result.

This audit measures actual exposure: how many tokens are removed per cell by
the treatment, and by each matched control set.

Direction of concern:
  treatment removes FEWER tokens than controls  -> null could be an artefact
  treatment removes AS MANY OR MORE             -> design biases toward
                                                   detecting harm, so a null
                                                   is conservative

Outputs (revision/outputs/):
  E9_token_occurrence.csv    -- per control set
  E9_token_occurrence.json   -- summary per treatment arm

Usage:  python E9_token_occurrence_audit.py [--dataset pbmc3k]
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

# Repository-relative paths. Scripts live in analysis/; everything they read
# and write is inside this repository.
REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUT = REPO / "outputs"
CACHE = REPO / "cache"
for _d in (OUT, CACHE):
    _d.mkdir(parents=True, exist_ok=True)
BASE = REPO  # legacy alias




def tokens_removed_per_cell(tokenised, token_ids):
    """Mean number of tokens a deletion set removes from each cell."""
    tset = set(int(t) for t in token_ids)
    return np.array([sum(1 for tok in cell if int(tok) in tset)
                     for cell in tokenised], dtype=float)


def audit(dataset="pbmc3k"):
    tok_path = CACHE / f"E2_{dataset}_tokenized.json"
    if not tok_path.exists():
        print(f"  No tokenised cache for {dataset}; run --setup first.")
        return None

    with open(tok_path) as f:
        tokenised = json.load(f)["tokenized"]

    treatment = pd.read_csv(OUT / "E2_treatment_genes.csv")
    arms = [
        ("full", treatment["token_id"].tolist(),
         CACHE / f"E2_matched_controls_{dataset}.json"),
        ("no_ribo_mito",
         treatment[~treatment["gene_class"].isin(
             ["ribosomal", "mitochondrial"])]["token_id"].tolist(),
         CACHE / f"E2_matched_controls_{dataset}_no_ribo_mito.json"),
    ]

    rows, summary = [], {}
    print("=" * 74)
    print(f"  E9 - Token-occurrence audit ({dataset})")
    print("=" * 74)
    print(f"  {len(tokenised)} cells\n")

    for arm, treat_tokens, ctrl_path in arms:
        if not ctrl_path.exists():
            print(f"  [{arm}] no matched controls; skipping")
            continue
        with open(ctrl_path) as f:
            controls = json.load(f)

        t_per_cell = tokens_removed_per_cell(tokenised, treat_tokens)
        t_mean = float(t_per_cell.mean())

        c_means = []
        for i, cset in enumerate(controls):
            m = float(tokens_removed_per_cell(tokenised, cset).mean())
            c_means.append(m)
            rows.append(dict(dataset=dataset, arm=arm, control_idx=i,
                             mean_tokens_removed=m))
        c_means = np.array(c_means)
        z = (t_mean - c_means.mean()) / c_means.std(ddof=1)
        frac_ge = float((c_means >= t_mean).mean())

        verdict = ("treatment at least as exposed as controls -> a null is "
                   "conservative" if z >= 0 else
                   "treatment LESS exposed than controls -> null may be an "
                   "exposure artefact")
        print(f"  [{arm}]  n_genes={len(treat_tokens)}  "
              f"n_controls={len(controls)}")
        print(f"      treatment removes {t_mean:.2f} tokens/cell")
        print(f"      controls  remove  {c_means.mean():.2f} "
              f"(SD {c_means.std(ddof=1):.2f})")
        print(f"      z = {z:+.2f}   {frac_ge * 100:.0f}% of control sets "
              f"remove at least as many")
        print(f"      -> {verdict}\n")

        summary[arm] = dict(
            n_treatment_genes=len(treat_tokens), n_controls=len(controls),
            treatment_tokens_per_cell=t_mean,
            control_tokens_per_cell_mean=float(c_means.mean()),
            control_tokens_per_cell_sd=float(c_means.std(ddof=1)),
            z=float(z), frac_controls_ge_treatment=frac_ge,
            interpretation=verdict)

    if rows:
        pd.DataFrame(rows).to_csv(OUT / "E9_token_occurrence.csv", index=False)
        with open(OUT / "E9_token_occurrence.json", "w") as f:
            json.dump({"dataset": dataset, "arms": summary,
                       "note": ("Exposure is measured on the tokenised cells "
                                "actually analysed, not on global expression "
                                "statistics. This is the quantity that "
                                "determines how much signal a deletion "
                                "removes.")}, f, indent=2)
        print(f"  Saved: {OUT/'E9_token_occurrence.csv'}")
        print(f"         {OUT/'E9_token_occurrence.json'}")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="pbmc3k")
    audit(ap.parse_args().dataset)
