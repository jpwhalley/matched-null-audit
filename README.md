# Matched-null auditing of gene-embedding outliers in single-cell foundation models

Reproduction package for the PSB 2027 paper of the same name (J. Whalley).

This repository exists to re-run the analyses reported in the paper. It holds
the code, the versioned inputs and their checksums, the saved outputs, and a
table mapping each reported result to the script and file that produced it.
It is not a general supplement and carries no findings or project history.

## Install

Requires Python 3.11 and [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/jpwhalley/matched-null-audit
cd matched-null-audit
uv sync
```

## Reproduce

```bash
./run_all.sh
```

About a minute. Needs only the shipped inputs: no model checkpoints, no GPU,
no network. It rebuilds all four figures and regenerates every analysis except
two, then checks the regenerated outputs against the committed ones and exits
non-zero on any difference.

The two it skips are the matched deletion test, which needs the Geneformer
checkpoint and about 18 hours, and the token-exposure audit, which reads the
tokenised cell matrices that are too large to ship. Their outputs are committed
so the reported numbers can be checked against them.

## Reproduce the expensive stages

Both need model checkpoints; see `DATA_MANIFEST.md` and the `D`-series
notebooks for acquisition.

```bash
RUN_ABLATION=1 ./run_all.sh        # matched deletion test and token audit, ~18 h
REGENERATE_ESM2=1 ./run_all.sh     # recompute ESM-2 embeddings, ~2 h
```

The deletion test checkpoints every ten control sets and resumes on restart.
The matched-control draws are shipped under `cache/`, so the null bands can be
recomputed without redrawing them; `E2_downstream_ablation.py` refuses any
cache whose specification sidecar disagrees with the settings in force.

Results depend on the compute device. The supervised probe baseline differs by
5.1e-4 macro-F1 between CPU and Metal Performance Shaders on an otherwise
identical stack, so each result records a full environment fingerprint and the
baseline is pinned to the device the ablation ran on.

## Layout

```
analysis/           analysis scripts, one per reported experiment
data/               shipped inputs and checksums
cache/              matched-control draws and their specification sidecars
outputs/            saved results
figures/            the four manuscript figures
notebooks/          D01-D04 data acquisition, P01 the geometry screen
run_all.sh          the pipeline above
DATA_MANIFEST.md    provenance, versions, licences, reproduction cost
MANUSCRIPT_TRACEABILITY.md   reported result -> script -> input -> output
```

## Citation

See `CITATION.cff`.

## Licence

MIT, see `LICENSE`. Model weights and external databases are not
redistributed here and carry their own terms; see `DATA_MANIFEST.md`.
