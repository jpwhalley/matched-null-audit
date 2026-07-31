"""
Clear the old matched-control caches so E2 --setup rebuilds them.

Refuses to run unless an archive manifest exists and verifies, because the
control caches are git-ignored and the archive is their only copy.

Deletes: E2_matched_controls_*  and  E2_ablation_ckpt_*
Keeps:   E2_baseline_*, E2_*_tokenized.json, E2_*_gene_stats.csv, *.h5ad

  python revision/notebooks/clear_controls.py
"""

import hashlib
import json
from pathlib import Path

_HERE = Path(__file__).resolve()
if _HERE.parent.name == "analysis":
    BASE = _HERE.parent.parent
    CACHE = BASE / "cache"
else:
    BASE = _HERE.parent.parent.parent
    CACHE = BASE / "revision" / "cache"

MANIFEST = BASE / "archive" / "ribo_v1_digit_anchored" / "MANIFEST.json"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    if not MANIFEST.exists():
        raise SystemExit(
            f"No archive manifest at {MANIFEST}\n"
            f"Run this first, it is the only copy of these controls:\n"
            f"  python revision/notebooks/archive_controls.py")
    man = json.loads(MANIFEST.read_text())
    bad = [e["archived_as"] for e in man["files"]
           if not (MANIFEST.parent / e["archived_as"]).exists()
           or sha256(MANIFEST.parent / e["archived_as"]) != e["sha256"]]
    if bad:
        raise SystemExit(f"Archive is incomplete or corrupt: {bad}\n"
                         f"Refusing to delete anything.")
    print(f"  archive verified: {len(man['files'])} files")

    targets = sorted(p for p in CACHE.iterdir() if p.is_file()
                     and (p.name.startswith("E2_matched_controls_")
                          or p.name.startswith("E2_ablation_ckpt_")))
    if not targets:
        print("  nothing to clear; cache is already rebuilt-ready")
    for p in targets:
        p.unlink()
        print(f"  deleted  {p.name}")

    left = [p.name for p in CACHE.iterdir() if p.is_file()
            and (p.name.startswith("E2_matched_controls_")
                 or p.name.startswith("E2_ablation_ckpt_"))]
    if left:
        raise SystemExit(f"STILL PRESENT after deletion: {left}\n"
                         f"Something is restoring them. Pause Google Drive "
                         f"syncing on this folder and run again.")
    print("  cache is clear")

    kept = sorted(p.name for p in CACHE.iterdir() if p.is_file()
                  and p.name.startswith(("E2_baseline_", "E2_pbmc3k_",
                                         "E2_tabula_sapiens_")))
    print(f"  kept {len(kept)} baseline/tokenisation files:")
    for k in kept:
        print(f"     {k}")


if __name__ == "__main__":
    main()
