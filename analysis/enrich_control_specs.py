"""
Enrich the matched-control specification sidecars.

WHY THIS EXISTS
---------------
The first version of the sidecar recorded only the ribosomal panel provenance
and the draw count. That is under-specified: the primary and the sensitivity
sidecars for a dataset came out byte-identical, so the guard in
`E2_downstream_ablation.py` could not tell one arm's cache from the other's.
The artefacts themselves were correct -- the cache-to-ablation chain was
verified directly, 50 genes and 100 evaluations for the Tabula Sapiens primary
against 36 genes and none for its sensitivity arm -- but a guard that cannot
distinguish the two is a guard waiting to pass something it should refuse.

This script adds the distinguishing fields. It reads the caches that already
exist and rebuilds nothing: no controls are drawn, no analysis is rerun, and
the recorded panel provenance is preserved untouched.

Added per sidecar:
  dataset, arm, n_draws, genes_per_draw,
  treatment_gene_list_sha256, control_cache_sha256

The treatment-gene-list hash is over the arm's sorted Ensembl IDs, so it is
stable against row order and against symbol renaming.

Usage:
  python enrich_control_specs.py            # report only
  python enrich_control_specs.py --apply
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
if _HERE.parent.name == "analysis":
    BASE = _HERE.parent.parent
    CACHE, OUT = BASE / "cache", BASE / "outputs"
else:
    BASE = _HERE.parent.parent.parent
    CACHE, OUT = BASE / "revision" / "cache", BASE / "revision" / "outputs"

TREATMENT = OUT / "E2_treatment_genes.csv"
EXCLUDED_FROM_SENSITIVITY = {"ribosomal", "mitochondrial"}

ARMS = [("primary", ""), ("sensitivity", "_no_ribo_mito")]


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def arm_gene_hash(df: pd.DataFrame, arm: str):
    d = df if arm == "primary" else df[~df.gene_class.isin(EXCLUDED_FROM_SENSITIVITY)]
    ids = sorted(d.ensembl_id.astype(str))
    return sha256_bytes("\n".join(ids).encode()), len(ids)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    if not TREATMENT.exists():
        raise SystemExit(f"Need {TREATMENT} to hash the treatment gene list.")
    tg = pd.read_csv(TREATMENT)

    datasets = sorted({p.name.split("E2_matched_controls_")[1]
                       .split("_no_ribo_mito")[0].removesuffix(".json")
                       for p in CACHE.glob("E2_matched_controls_*.json")
                       if not p.name.endswith(("_spec.json", "_detail.json"))})
    if not datasets:
        raise SystemExit(f"No control caches under {CACHE}")

    changed = 0
    for ds in datasets:
        for arm, suffix in ARMS:
            cache_p = CACHE / f"E2_matched_controls_{ds}{suffix}.json"
            spec_p = CACHE / f"E2_matched_controls_{ds}{suffix}_spec.json"
            if not cache_p.exists() or not spec_p.exists():
                print(f"  {ds:16s} {arm:12s} absent, skipped")
                continue

            raw = cache_p.read_bytes()
            draws = json.loads(raw)
            sizes = {len(x) for x in draws}
            gh, n_expected = arm_gene_hash(tg, arm)

            spec = json.loads(spec_p.read_text())
            add = {
                "dataset": ds,
                "arm": arm,
                "n_draws": len(draws),
                "genes_per_draw": sorted(sizes)[0] if len(sizes) == 1 else sorted(sizes),
                "treatment_gene_list_sha256": gh,
                "control_cache_sha256": sha256_bytes(raw),
            }
            if len(sizes) != 1:
                print(f"  {ds:16s} {arm:12s} WARNING ragged draw sizes {sizes}")
            elif sorted(sizes)[0] != n_expected:
                print(f"  {ds:16s} {arm:12s} WARNING draw size "
                      f"{sorted(sizes)[0]} != treatment list {n_expected}")

            merged = {**spec, **add}
            same = merged == spec
            print(f"  {ds:16s} {arm:12s} n={add['n_draws']:4d} "
                  f"genes/draw={add['genes_per_draw']}  "
                  f"tx={gh[:12]}  cache={add['control_cache_sha256'][:12]}"
                  f"{'  (already current)' if same else ''}")
            if args.apply and not same:
                spec_p.write_text(json.dumps(merged, indent=2) + "\n")
                changed += 1

    # The point of the exercise: no two sidecars may now be identical.
    specs = sorted(CACHE.glob("E2_matched_controls_*_spec.json"))
    by_hash = {}
    for p in specs:
        by_hash.setdefault(sha256_bytes(p.read_bytes()), []).append(p.name)
    dupes = {h: v for h, v in by_hash.items() if len(v) > 1}
    print()
    if dupes:
        for h, v in dupes.items():
            print(f"  STILL IDENTICAL {h[:12]}: {v}")
        if args.apply:
            sys.exit(1)
    else:
        print(f"  all {len(specs)} sidecars are now distinct")
    print("  (report only; rerun with --apply to write)" if not args.apply
          else f"  wrote {changed} sidecar(s)")


if __name__ == "__main__":
    main()
