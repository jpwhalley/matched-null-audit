#!/usr/bin/env bash
# Reproduce the analyses reported in the manuscript.
#
# Fails loudly on any missing script or non-zero exit. Does not silently skip.
#
# The default path needs only the shipped inputs: no model checkpoints, no GPU,
# no network. It takes about a minute. Two stages are expensive and opt-in:
#
#   REGENERATE_ESM2=1   recompute ESM-2 embeddings from the checkpoint (~2 h)
#   RUN_ABLATION=1      run the matched deletion test (~18 h)
#
# Data acquisition for those paths is in notebooks/D01-D04; see DATA_MANIFEST.md.

set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-uv run python}"

run () {
  local script="analysis/$1"; shift
  [[ -f "$script" ]] || { echo "MISSING: $script" >&2; exit 1; }
  echo; echo "=== $script $* ==="
  $PY "$script" "$@"
}

echo "### Cross-model agreement — Table 1"
run E10_cross_model_agreement.py

echo "### Caller robustness — Figure 2"
run E3_outlier_robustness.py

echo "### Protein-sequence control — Figure 3"
# The ESM-2 geometry table is shipped, so the default verifies it against
# data/CHECKSUMS.json and recomputes the comparison from it.
if [[ "${REGENERATE_ESM2:-0}" == "1" ]]; then
  run E1_stage1_mapping.py --all
  run E1_stage2_esm2.py --all
else
  run E1_stage2_esm2.py --verify-shipped
fi

echo "### Covariate-adjusted association — Table 2"
run E8_clinvar_adjusted.py
run E6_class_association.py

echo "### Exploratory dependency analysis"
run E11_dependency.py

echo "### Matched deletion test — Figure 4  (expensive)"
if [[ "${RUN_ABLATION:-0}" != "1" ]]; then
  echo "  Skipped. Set RUN_ABLATION=1 to run the ablation and token audit."
else
  run E2_downstream_ablation.py --setup --baseline --ablation --evaluate --datasets pbmc3k
  # The matched-control draws are shipped under cache/, but the tokenised cells
  # they index (E2_<dataset>_tokenized.json) are not, so the token audit can
  # only run after --setup has rebuilt them.
  run E9_token_occurrence_audit.py --dataset pbmc3k
fi

# Reads only shipped outputs and takes seconds, so it runs on every invocation
# and stays under the drift gate below.
echo "### Clustering-metric diagnostic"
run E7_cluster_metric_diagnostic.py

echo "### Figures"
run make_psb_figures.py

echo; echo "Done. Outputs in outputs/, figures in figures/."
echo
echo "Verifying nothing drifted:"
if command -v git >/dev/null && git rev-parse --git-dir >/dev/null 2>&1; then
  # outputs/ is gated strictly: floats are serialised at a documented precision
  # (analysis/_precision.py) so library-version noise cannot trip it.
  # figures/ is advisory: PDF bytes depend on the rendering stack.
  ok=0
  if git diff --exit-code --stat -- outputs; then
    echo "  outputs: clean."
  else
    echo "  DRIFT in outputs/ — regenerated results differ from the commit." >&2
    ok=1
  fi
  if ! git diff --quiet -- figures; then
    echo "  NOTE: figures/ differ byte-wise. Expected on a different rendering"
    echo "        stack; check visually rather than by hash."
  fi
  exit $ok
fi
