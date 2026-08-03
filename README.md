# Do Geometric Outliers Identify Important Genes in Single-Cell Foundation Models?

Reproduction package for the PSB 2027 paper of the same name (J. P. Whalley).

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
no network. Starting from the checksum-pinned geometry tables, the ESM-2
geometry table and the saved deletion outputs committed here, it regenerates
the downstream analyses and all five manuscript figures, then checks the
regenerated outputs against the committed ones and exits non-zero on any
difference.

It is a reproduction from frozen analysis inputs, not a rebuild from raw
checkpoints. By default it does not re-extract the gene-embedding matrices
from the models, recompute the ESM-2 embeddings, or rerun the matched
deletion test and the token-exposure audit. The deletion test needs the
Geneformer checkpoint and about 18 hours; the token-exposure audit reads
tokenised cell matrices too large to ship. The outputs of all of these are
committed, so the reported numbers can be checked against them, and each
stage can be re-run from source with the opt-in commands below.

## Reproduce the expensive stages

Both need model checkpoints; see `DATA_MANIFEST.md` and the `D`-series
notebooks for acquisition.

```bash
RUN_ABLATION=1 ./run_all.sh        # PBMC3k matched deletion and token audit, ~18 h
REGENERATE_ESM2=1 ./run_all.sh     # recompute ESM-2 embeddings, ~2 h
```

`RUN_ABLATION=1` covers PBMC3k only. The Tabula Sapiens replication is a
separate run with its own settings and is not reached by any environment
variable:

```bash
uv run python analysis/E2_downstream_ablation.py \
    --setup --baseline --ablation --evaluate \
    --datasets tabula_sapiens --primary-only --subsample 4000 --n-bootstrap 100
uv run python analysis/E9_token_occurrence_audit.py --dataset tabula_sapiens
```

`--primary-only` runs the top-*k* treatment and its matched null and skips the
sensitivity arm and the *k*-sweep, which is what the paper reports for this
dataset.

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
figures/            the five manuscript figures
notebooks/          D01, D02 and D04 data acquisition; P01 the geometry screen
run_all.sh          the pipeline above
DATA_MANIFEST.md    provenance, versions, licences, reproduction cost
MANUSCRIPT_TRACEABILITY.md   reported result -> script -> input -> output
```

## Citation

See `CITATION.cff`.

## Licence

MIT, see `LICENSE`. Model weights and external databases are not
redistributed here and carry their own terms; see `DATA_MANIFEST.md`.
