#!/usr/bin/env bash
# Reproduce the analyses reported in the manuscript.
#
# Fails loudly on any missing script or non-zero exit. Does NOT silently skip.
# Data acquisition (notebooks/D01-D04) must be run first; see DATA_MANIFEST.md.

set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-uv run python}"

run () {
  local script="analysis/$1"; shift
  [[ -f "$script" ]] || { echo "MISSING: $script" >&2; exit 1; }
  echo; echo "=== $script $* ==="
  $PY "$script" "$@"
}

echo "### Stage 2 — caller robustness"
run E3_outlier_robustness.py

echo "### Stage 3 — non-transcriptomic control"
# E1 stage 2 regenerates ESM-2 embeddings from the checkpoint (~2 h, GPU
# advisable). The resulting geometry table is SHIPPED, so by default we verify
# against it rather than recompute. Set REGENERATE_ESM2=1 to recompute.
if [[ "${REGENERATE_ESM2:-0}" == "1" ]]; then
  run E1_stage1_mapping.py --all
  run E1_stage2_esm2.py --all
else
  echo "  Using shipped outputs/E1_esm2_geometry.csv."
  echo "  Set REGENERATE_ESM2=1 to recompute from the ESM-2 checkpoint."
  [[ -f outputs/E1_esm2_geometry.csv ]] || { echo "MISSING: outputs/E1_esm2_geometry.csv" >&2; exit 1; }
fi

echo "### Stage 5 — covariate-aware annotation"
run E8_clinvar_adjusted.py
run E6_class_association.py

echo "### Stage 4 — matched deletion test  (SLOW: ~18 h, checkpoints every 10 controls)"
if [[ "${SKIP_ABLATION:-0}" == "1" ]]; then
  echo "  SKIP_ABLATION=1 set; skipping."
else
  run E2_downstream_ablation.py --setup --baseline --ablation --evaluate --datasets pbmc3k
  run E9_token_occurrence_audit.py --dataset pbmc3k
  run E7_cluster_metric_diagnostic.py
fi

echo "### Figures"
run make_psb_figures.py

echo; echo "Done. Outputs in outputs/, figures in figures/."
echo
echo "Verifying nothing drifted:"
if command -v git >/dev/null && git rev-parse --git-dir >/dev/null 2>&1; then
  if git diff --exit-code --stat -- outputs figures; then
    echo "  clean: regenerated outputs match the committed ones."
  else
    echo "  DRIFT: regenerated outputs differ from the commit (see above)." >&2
    exit 1
  fi
fi
