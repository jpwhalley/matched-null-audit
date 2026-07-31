"""
Archive the digit-anchored (pre-2026-07-30) matched controls and their results.

WHY THIS EXISTS
---------------
The ribosomal panel was replaced on 2026-07-30 with a pinned panel built from
HGNC gene groups 728, 729 and 646. The prior symbol regex both omitted curated
members and admitted non-ribosomal genes, including RPS6K signalling kinases.
Because those labels affect the matched-control strata, the eligible pool
changes and the matched nulls must be rebuilt.

The original controls are not a mistake to be deleted. Under pre-commitment 4 in
`prespecification/results_addendum_2026-07-29.md` they are retained as a
disclosed provenance and sensitivity analysis. The control caches are
git-ignored, so without this step a wildcard `rm` would destroy the only copy.

This script COPIES, hashes and manifests. It deletes nothing; the removal
command is printed for you to run deliberately.

Usage:
  python archive_controls.py            # archive, then print the rm command
  python archive_controls.py --verify   # re-check an existing archive
"""

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

_HERE = Path(__file__).resolve()
if _HERE.parent.name == "analysis":
    BASE = _HERE.parent.parent
    OUT, CACHE = BASE / "outputs", BASE / "cache"
else:
    BASE = _HERE.parent.parent.parent
    OUT, CACHE = BASE / "revision" / "outputs", BASE / "revision" / "cache"

LABEL = "ribo_v1_digit_anchored"
ARCHIVE = BASE / "archive" / LABEL

PATTERNS = [
    (CACHE, "E2_matched_controls_*.json"),
    (CACHE, "E2_ablation_ckpt_*.json"),
    (OUT,   "E2_ablation_*.json"),
    (OUT,   "E2_verdict_*.json"),
    (OUT,   "E2_matched_controls_*_balance.csv"),
    (OUT,   "E9_token_occurrence*.json"),
    (OUT,   "E9_token_occurrence*.csv"),
    (OUT,   "E7_cluster_metric_diagnostic.*"),
]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def collect():
    seen, files = set(), []
    for root, pat in PATTERNS:
        if not root.exists():
            continue
        for f in sorted(root.glob(pat)):
            if f.is_file() and f not in seen:
                seen.add(f)
                files.append(f)
    return files


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()
    manifest_path = ARCHIVE / "MANIFEST.json"

    if args.verify:
        if not manifest_path.exists():
            raise SystemExit(f"No manifest at {manifest_path}")
        man = json.loads(manifest_path.read_text())
        bad = 0
        for e in man["files"]:
            a = ARCHIVE / e["archived_as"]
            if not a.exists():
                print(f"  MISSING  {e['archived_as']}"); bad += 1
            elif sha256(a) != e["sha256"]:
                print(f"  CORRUPT  {e['archived_as']}"); bad += 1
        print(f"  {len(man['files']) - bad}/{len(man['files'])} verified")
        raise SystemExit(1 if bad else 0)

    files = collect()
    if not files:
        raise SystemExit(
            "Nothing to archive. The caches are git-ignored and empty in a "
            "clean clone; run this on the compute workspace.")

    ARCHIVE.mkdir(parents=True, exist_ok=True)
    entries = []
    print(f"  Archiving {len(files)} files to {ARCHIVE}")
    for f in files:
        rel = f.relative_to(BASE)
        dest_name = str(rel).replace("/", "__")
        shutil.copy2(f, ARCHIVE / dest_name)
        digest = sha256(f)
        entries.append(dict(
            source=str(rel), archived_as=dest_name, sha256=digest,
            bytes=f.stat().st_size,
            mtime=datetime.fromtimestamp(
                f.stat().st_mtime, timezone.utc).isoformat()))
        print(f"    {digest[:16]}  {rel}")

    manifest_path.write_text(json.dumps({
        "label": LABEL,
        "reason": (
            "Matched controls and deletion results produced under the "
            "digit-anchored ribosomal pattern ^(RPL|RPS|MRPL|MRPS)\\d, which "
            "omitted curated ribosomal genes and admitted non-ribosomal "
            "RPS6K kinases. Retained as the disclosed provenance and "
            "sensitivity analysis under pre-commitment 4 of "
            "prespecification/results_addendum_2026-07-29.md."),
        "archived_utc": datetime.now(timezone.utc).isoformat(),
        "n_files": len(entries),
        "files": entries,
    }, indent=2))

    print(f"\n  Manifest: {manifest_path}")
    print("  Verify with:  python archive_controls.py --verify")
    print("\n  Once verified, remove ONLY these to force a clean rebuild.")
    print("  Use find, not rm with globs: in zsh an unmatched glob aborts the")
    print("  whole command line, so a missing checkpoint pattern would silently")
    print("  leave the control caches in place.")
    print(f"    find {CACHE} -maxdepth 1 -name 'E2_matched_controls_*.json' -delete")
    print(f"    find {CACHE} -maxdepth 1 -name 'E2_ablation_ckpt_*.json' -delete")
    print("  Keep E2_baseline_*, E2_*_tokenized.json and E2_*_gene_stats.csv:")
    print("  baselines and treatment selection are not being recomputed.")


if __name__ == "__main__":
    main()
