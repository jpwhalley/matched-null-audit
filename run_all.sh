#!/usr/bin/env bash
# Reproduce the analyses reported in the manuscript.
#
# Fails loudly on any missing script or non-zero exit. Does NOT silently skip.
#
# The default path needs only the SHIPPED inputs -- no model checkpoints, no
# GPU, no network. Data acquisition (notebooks/D01-D04) is required only for
# full regeneration: REGENERATE_ESM2=1 and the E2 ablation. See DATA_MANIFEST.md.

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
  # Verifies the shipped geometry against data/CHECKSUMS.json, then recomputes
  # the scFM comparison and verdict from it. ~1.5 s.
  run E1_stage2_esm2.py --verify-shipped
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
  # outputs/ is gated strictly: floats are serialised at a documented
  # precision (analysis/_precision.py) so library-version noise cannot trip it.
  # figures/ is advisory: PDF bytes depend on the font and rendering stack even
  # with timestamps suppressed, so a difference there is reported, not fatal.
  ok=0
  if git diff --exit-code --stat -- outputs; then
    echo "  outputs: clean."
  else
    echo "  DRIFT in outputs/ — regenerated results differ from the commit." >&2
    ok=1
  fi
  if ! git diff --quiet -- figures; then
    echo "  NOTE: figures/ differ byte-wise. Expected on a different font or"
    echo "        rendering stack; check visually rather than by hash."
  fi
  exit $ok
fi
