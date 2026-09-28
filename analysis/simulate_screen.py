"""
Screen-calibration simulation for the PSB 2027 camera-ready (reviewer 3, R3.3).

Question answered: how does the four-metric outlier screen of Section 2.1
behave on matrices with known structure? Two readouts only: (1) the call rate
under nulls with no planted structure, per metric and for the union rule;
(2) recovery of planted outliers, reported as a firing matrix (which metrics
call each planted class at each perturbation size), on matrices where all
screen statistics are recomputed exactly as the real screen computes them.
It simulates the screen, not biology, and says nothing about the three real
models beyond calibrating the instrument applied to them.

Design (approved by Justin, 2026-09-10, with his four clarifications):
  * n = 20,000 genes; d fixed at 256, 512 and 768; seeds fixed by a timing
    benchmark of one 20,000 x 768 replicate (the benchmark sets only the seed
    count, never the dimensions).
  * Null A: isotropic Gaussian. Null B: Null A plus a shared mean direction
    u with magnitude a = c * sqrt(d / (1 - c^2)), so the mean cosine to the
    centroid is c = 0.243, the median across the three models of each model's
    mean cosine to its centroid (Geneformer 0.290, scGPT 0.243,
    scFoundation 0.009; repository data/*gene_embedding_geometry.csv,
    special tokens excluded).
  * Planted outliers are planted on Null B only, one mechanism per planted
    matrix (so that one mechanism's extremes cannot inflate the z-score
    denominators of another's), 100 genes per matrix, three perturbation
    sizes chosen to bracket the |z| > 3 threshold (about 2, 4 and 8 null
    standard deviations of the targeted metric): norm (x -> s x; s = 1.05,
    1.10, 1.20) and direction (rotate x away from u in the plane of u and
    x's orthogonal component by 5, 10, 20 degrees, norm preserved).
    Isolation is not planted: the attempted construction (a norm-preserving
    jitter orthogonal to the shared direction; pilot v1, 2026-09-10, run with
    the Euclidean isolation metric of that pilot) did not produce an
    isolation-only perturbation, moving the cosine and leaving the isolation
    z-score unchanged, so the isolation metric is instead characterised by
    its null call rate and by its co-firing on the planted genes (under the
    cosine metric it co-fires on direction-planted genes and, being
    scale-invariant, never on norm-planted genes). The same 100 planted genes are used for every planted matrix of a
    seed (paired recovery curves across sizes and mechanisms). The firing
    matrix (which metrics call the planted genes) and the call rate among
    unplanted genes on each planted matrix (masking) are the outputs.
  * Reporting: results are summarised separately by dimension, as seed-level
    means and ranges across seeds, never by pooling genes across seeds or
    dimensions. The fixed factors and angles are different standardised
    effects at different d (the "about 2, 4 and 8 null standard deviations"
    description holds for the d = 768 pilot only); the mean z per metric is
    written for every planted group so the standardised effect at each d is
    visible.
  * Seeds: seed 0 was the timing and design pilot and is excluded from the
    full run, which uses seeds 1 to 20.
  * Screen: embedding norm, Euclidean distance to the vocabulary centroid,
    cosine to the centroid, mean cosine distance to the 10 nearest
    neighbours (self excluded), exactly as notebooks/P01_embedding_geometry
    computes them (NearestNeighbors(metric='cosine')); each z-scored across
    all genes of the matrix
    being screened (planted genes included, as the real screen would);
    outlier if |z| > 3 on any metric. Gaussian reference points (0.27% per
    metric, at most ~1.1% for the union) are references, not expectations.

Usage:
  python simulate_screen.py --benchmark                      # pilot: seed 0, d=768, timed
  python simulate_screen.py --dims 768 --seed-from 1 --seed-to 20   # part of the full run
  python simulate_screen.py --summarise                      # merge parts, write summary
Version note (2026-09-16, camera-ready v2): the 2026-09-10 run used mean
Euclidean distance to the 10 nearest neighbours for isolation, which is not
the metric of the real screen; the function was corrected to cosine distance
and seeds 1-20 were rerun with the design otherwise unchanged. Null union
rates moved from 0.6-0.7% to 0.8-0.9% (isolation null calls 0.27-0.35%);
recovery changed by at most about two percentage points at any size.

Full-run parts go to <out>/simulate_screen_runs/d<d>_s<seed>.csv (restart-safe);
--summarise merges them into simulate_screen_run_results.csv, writes
simulate_screen_run_summary.csv (seed-level mean, min, max by dimension,
matrix, group and metric) and simulate_screen_run_config.json. Exact nearest
neighbours by chunked matrix products; numpy only.
"""
import argparse, json, time, sys
from pathlib import Path
import numpy as np

N_GENES = 20_000
DIMS = (256, 512, 768)
C_COS = 0.243
N_PLANT = 100
SIZES = {"norm": (1.05, 1.10, 1.20), "direction": (5.0, 10.0, 20.0)}
K = 10
Z_CUT = 3.0
METRICS = ("norm", "dist_centroid", "cos_centroid", "isolation")


def knn_mean_dist(X, k=K, chunk=1000):
    """Mean cosine distance to the k nearest neighbours (cosine metric), self excluded. Exact.
    Matches the real screen (P01_embedding_geometry.ipynb: NearestNeighbors(metric='cosine'))."""
    n = X.shape[0]
    Xn = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    out = np.empty(n)
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        D = 1.0 - Xn[s:e] @ Xn.T
        D[np.arange(e - s), np.arange(s, e)] = np.inf  # exclude self
        part = np.partition(D, k, axis=1)[:, :k]
        out[s:e] = part.mean(axis=1)
    return out


def screen(X):
    """The four metrics, z-scored across all rows; returns z (n x 4) and calls."""
    centroid = X.mean(axis=0)
    norm = np.linalg.norm(X, axis=1)
    diff = X - centroid
    dist = np.linalg.norm(diff, axis=1)
    cos = (X @ centroid) / (norm * np.linalg.norm(centroid) + 1e-12)
    iso = knn_mean_dist(X)
    M = np.stack([norm, dist, cos, iso], axis=1)
    z = (M - M.mean(axis=0)) / M.std(axis=0, ddof=0)
    calls = np.abs(z) > Z_CUT
    return z, calls


def null_a(rng, n, d):
    return rng.standard_normal((n, d))


def shared_direction(rng, d):
    u = rng.standard_normal(d); return u / np.linalg.norm(u)


def null_b(rng, n, d, c=C_COS):
    a = c * np.sqrt(d / (1.0 - c * c))
    u = shared_direction(rng, d)
    return null_a(rng, n, d) + a * u, u, a


def plant(rng, X, u, mechanism, level, idx):
    """Return a planted copy of X for one mechanism and size; idx is the shared planted set."""
    X = X.copy()
    if mechanism == "norm":
        X[idx] *= SIZES["norm"][level]
    elif mechanism == "direction":
        theta = np.deg2rad(SIZES["direction"][level])
        for i in idx:
            x = X[i]; nx = np.linalg.norm(x)
            par = x @ u; perp = x - par * u; nperp = np.linalg.norm(perp)
            v = perp / nperp
            new = np.arctan2(nperp, par) + theta  # angle to u, increased by theta
            X[i] = nx * (np.cos(new) * u + np.sin(new) * v)
    else:
        raise ValueError(mechanism)
    return X, idx


def summarise(rows, tag, d, seed, matrix, z, calls, idx=None, mechanism=None, level=None, size=None):
    n = z.shape[0]
    union = calls.any(axis=1)
    planted = np.zeros(n, bool)
    if idx is not None: planted[idx] = True
    groups = {"all": np.ones(n, bool)}
    if idx is not None:
        groups["unplanted"] = ~planted; groups["planted"] = planted
    for g, mask in groups.items():
        rec = {"tag": tag, "d": d, "seed": seed, "matrix": matrix, "mechanism": mechanism, "level": level,
               "size": size, "group": g, "n": int(mask.sum())}
        for j, met in enumerate(METRICS):
            rec[f"call_{met}"] = float(calls[mask, j].mean())
            rec[f"mean_z_{met}"] = float(z[mask, j].mean())
        rec["call_union"] = float(union[mask].mean())
        rows.append(rec)


def run(dims, seeds, out_dir, tag):
    rows = []; timings = []
    for d in dims:
        for seed in seeds:
            rng = np.random.default_rng(10_000 * d + seed)
            t0 = time.time()
            XA = null_a(rng, N_GENES, d); z, c = screen(XA); summarise(rows, tag, d, seed, "nullA", z, c)
            XB, u, a = null_b(rng, N_GENES, d); z, c = screen(XB); summarise(rows, tag, d, seed, "nullB", z, c)
            realised_cos = float(((XB @ u) / np.linalg.norm(XB, axis=1)).mean())
            idx = rng.choice(N_GENES, N_PLANT, replace=False)  # one planted set per seed: paired curves
            for mech in SIZES:
                for level in range(3):
                    XP, idx = plant(rng, XB, u, mech, level, idx)
                    z, c = screen(XP)
                    summarise(rows, tag, d, seed, f"{mech}_L{level + 1}", z, c, idx, mech, level + 1, SIZES[mech][level])
            dt = time.time() - t0
            timings.append({"d": d, "seed": seed, "seconds": round(dt, 1), "a": round(float(a), 3), "realised_mean_cos_nullB": round(realised_cos, 4)})
            print(f"d={d} seed={seed}: {dt:.1f} s (8 screens); a={a:.2f}, realised mean cos on Null B = {realised_cos:.3f}", flush=True)
    return rows, timings


def write_csv(path, rows):
    import csv
    keys = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n"); w.writeheader(); w.writerows(rows)


def summarise_runs(out_dir):
    """Merge per-(d, seed) parts; seed-level mean/min/max by dimension, matrix, group, metric."""
    import csv, glob
    parts = sorted(glob.glob(str(out_dir / "simulate_screen_runs" / "d*_s*.csv")))
    rows = []
    for pth in parts:
        with open(pth) as f:
            rows.extend(csv.DictReader(f))
    if not rows:
        raise SystemExit("no run parts found")
    write_csv(out_dir / "simulate_screen_run_results.csv", rows)
    val_cols = [k for k in rows[0] if k.startswith("call_") or k.startswith("mean_z_")]
    groups = {}
    for r in rows:
        key = (int(r["d"]), r["matrix"], r["group"])
        groups.setdefault(key, []).append(r)
    summary = []
    for (d, matrix, group), rs in sorted(groups.items()):
        rec = {"d": d, "matrix": matrix, "group": group, "n_seeds": len(rs),
               "seeds": ",".join(sorted({r["seed"] for r in rs}, key=int))}
        for c in val_cols:
            vals = np.array([float(r[c]) for r in rs])
            rec[f"{c}_mean"] = float(vals.mean()); rec[f"{c}_min"] = float(vals.min()); rec[f"{c}_max"] = float(vals.max())
        summary.append(rec)
    write_csv(out_dir / "simulate_screen_run_summary.csv", summary)
    seeds = sorted({int(r["seed"]) for r in rows}); dims = sorted({int(r["d"]) for r in rows})
    cfg = {"n_genes": N_GENES, "dims": dims, "seeds": seeds, "pilot_seed_excluded": 0, "c_cos": C_COS,
           "n_plant_per_mechanism": N_PLANT, "paired_planted_set_per_seed": True, "sizes": SIZES, "k": K,
           "z_cut": Z_CUT, "metrics": METRICS, "isolation_metric": "mean cosine distance to k nearest neighbours",
           "n_parts": len(parts), "numpy": np.__version__,
           "python": sys.version.split()[0]}
    with open(out_dir / "simulate_screen_run_config.json", "w") as f:
        json.dump(cfg, f, indent=2)
    print(f"merged {len(parts)} parts, {len(rows)} rows; dims {dims}; seeds {seeds[0]}..{seeds[-1]} (n={len(seeds)})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", action="store_true", help="pilot: one d=768 replicate with seed 0, timed")
    ap.add_argument("--dims", type=int, nargs="*", default=list(DIMS))
    ap.add_argument("--seed-from", type=int, default=1)
    ap.add_argument("--seed-to", type=int, default=20)
    ap.add_argument("--summarise", action="store_true")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "outputs")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.summarise:
        summarise_runs(args.out); return
    if args.benchmark:
        rows, timings = run([768], [0], args.out, "benchmark")
        write_csv(args.out / "simulate_screen_benchmark_results.csv", rows)
        with open(args.out / "simulate_screen_benchmark_config.json", "w") as f:
            json.dump({"pilot_seed": 0, "dims": [768], "timings": timings, "sizes": SIZES}, f, indent=2)
        return
    if args.seed_from < 1:
        raise SystemExit("seed 0 is the pilot seed and is excluded from the full run; use --seed-from 1 or higher")
    runs = args.out / "simulate_screen_runs"; runs.mkdir(exist_ok=True)
    for d in args.dims:
        for seed in range(args.seed_from, args.seed_to + 1):
            part = runs / f"d{d}_s{seed:02d}.csv"
            if part.exists():
                print(f"skip {part.name} (exists)"); continue
            rows, timings = run([d], [seed], args.out, "run")
            write_csv(part, rows)
    print("parts written to", runs)


if __name__ == "__main__":
    main()
