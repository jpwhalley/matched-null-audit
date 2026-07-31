"""
E10 — Direct cross-model agreement between geometric outlier sets.

WHY THIS EXISTS
---------------
The manuscript's central claim is that the three models produce model-specific
outlier profiles, and that a gene which is geometrically exceptional in one
model is largely not exceptional in another. Until now that claim rested on
indirect evidence: different outlier counts, different class enrichments, and
different behaviour under re-calling. None of those is the direct measurement.

This script supplies the direct measurement. On the shared vocabulary in
Table_S1.csv it reports, for each model pair:

  * the observed overlap of the two outlier sets
  * the overlap expected under independence
  * the Jaccard index
  * the Spearman correlation between the two continuous anomaly scores
  * a one-sided hypergeometric p for the observed overlap

The set-level and score-level views answer different questions and are both
reported: Jaccard asks whether the called sets coincide, the correlation asks
whether the underlying orderings agree genome-wide. A pair can overlap far
above chance and still have a low Jaccard, which is exactly what Geneformer
and scGPT do.

Comparisons are restricted to the shared vocabulary. The three models have
vocabularies of 20,271, 60,694 and 19,264 genes, so a raw comparison would
confound genuine disagreement with vocabulary coverage.

Outputs (outputs/):
  E10_cross_model_agreement.csv   -- one row per model pair
  E10_cross_model_agreement.json  -- headline figures and provenance

Usage:
  python E10_cross_model_agreement.py
"""

import hashlib
import itertools
import json
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _precision import clean  # documented serialisation precision

_HERE = Path(__file__).resolve()
if _HERE.parent.name == "analysis":
    REPO = _HERE.parent.parent
    DATA, OUT = REPO / "data", REPO / "outputs"
else:
    REPO = _HERE.parent.parent.parent
    DATA, OUT = REPO / "notebooks" / "data", REPO / "revision" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

TABLE = DATA / "Table_S1.csv"

# Membership is read from the committed outlier_class column so that this
# script cannot drift from the sets used everywhere else in the paper.
MEMBERSHIP = {
    "Geneformer":   ["GF-only", "GF∩scGPT", "GF∩SF", "All three"],
    "scGPT":        ["scGPT-only", "GF∩scGPT", "All three"],
    "scFoundation": ["SF-only", "GF∩SF", "All three"],
}
SCORE = {
    "Geneformer":   "gf_anomaly_score",
    "scGPT":        "scgpt_anomaly_score",
    "scFoundation": "sf_anomaly_score",
}


def main() -> None:
    sha = hashlib.sha256(TABLE.read_bytes()).hexdigest()
    t = pd.read_csv(TABLE, low_memory=False)
    n = len(t)

    print("=" * 74)
    print("  E10 - Direct cross-model agreement")
    print("=" * 74)
    print(f"  Table: {TABLE.name}  sha256[:16]={sha[:16]}")
    print(f"  Shared vocabulary: {n} genes")

    sets = {m: t.outlier_class.isin(c) for m, c in MEMBERSHIP.items()}
    for m, v in sets.items():
        print(f"    {m:13s} {int(v.sum()):5d} outliers "
              f"({100 * v.mean():.2f}% of shared vocabulary)")

    rows = []
    print()
    print(f"  {'pair':30s}{'obs':>6s}{'exp':>8s}{'Jaccard':>9s}"
          f"{'rho':>8s}{'hyperg p':>12s}")
    print("  " + "-" * 71)
    for a, b in itertools.combinations(MEMBERSHIP, 2):
        inter = int((sets[a] & sets[b]).sum())
        na, nb = int(sets[a].sum()), int(sets[b].sum())
        union = na + nb - inter
        expected = na * nb / n
        ok = t[SCORE[a]].notna() & t[SCORE[b]].notna()
        rho, rho_p = stats.spearmanr(t.loc[ok, SCORE[a]], t.loc[ok, SCORE[b]])
        # one-sided: P(X >= inter) under the hypergeometric null
        hyp_p = float(stats.hypergeom.sf(inter - 1, n, na, nb))
        rows.append(dict(
            model_a=a, model_b=b, n_shared_vocabulary=n,
            n_outliers_a=na, n_outliers_b=nb,
            observed_overlap=inter, expected_overlap=expected,
            fold_enrichment=inter / expected if expected else float("nan"),
            jaccard=inter / union if union else float("nan"),
            spearman_rho=float(rho), spearman_p=float(rho_p),
            hypergeometric_p=hyp_p, n_scored_both=int(ok.sum()),
        ))
        print(f"  {a + ' vs ' + b:30s}{inter:6d}{expected:8.1f}"
              f"{inter / union:9.3f}{rho:8.3f}{hyp_p:12.2e}")

    three = int((sets["Geneformer"] & sets["scGPT"]
                 & sets["scFoundation"]).sum())
    union_all = int((sets["Geneformer"] | sets["scGPT"]
                     | sets["scFoundation"]).sum())
    print()
    print(f"  Three-way intersection: {three}")
    print(f"  Union of all outlier sets: {union_all}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "E10_cross_model_agreement.csv", index=False)

    payload = {
        "headline": (
            "The two models that represent genes as discrete tokens "
            "(Geneformer, scGPT) overlap far above chance but still agree "
            "only weakly in absolute terms; scGPT and scFoundation overlap "
            "at chance. Agreement is graded, not uniform."
        ),
        "shared_vocabulary": n,
        "outlier_counts": {m: int(v.sum()) for m, v in sets.items()},
        "pairs": clean(rows),
        "three_way_intersection": three,
        "union_all_outliers": union_all,
        "note": (
            "Jaccard and the score correlation answer different questions "
            "and can diverge: a pair may overlap far above chance and still "
            "have a low Jaccard. Both are reported."
        ),
        "provenance": {"table": str(TABLE.relative_to(REPO)),
                       "sha256_16": sha[:16]},
    }
    (OUT / "E10_cross_model_agreement.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False))

    print()
    print(f"  Saved: {OUT / 'E10_cross_model_agreement.csv'}")
    print(f"         {OUT / 'E10_cross_model_agreement.json'}")


if __name__ == "__main__":
    main()
