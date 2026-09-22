#!/usr/bin/env bash
# GPU smoke recipe — one baseline, few generations, no Vast.ai required.
# Science path: evolution_orchestrator.py on a single GPU host.
#
# Prerequisites:
#   - cascade conda env activated (or venv with requirements.txt)
#   - EVAL_CMD = cattle-prod or protenix on PATH
#   - PXDESIGN_CMD set if Gen≥1 generation is enabled
#   - Phase 1 structures / metadata for at least one confirmed baseline
#
# Example (A100-class GPU):
#   export CASCADE_WORKERS=1
#   export CASCADE_MAX_GENERATIONS=2
#   export CASCADE_VARIANTS_PER_GEN=2
#   export CASCADE_BASELINE_IDS=NZ_JAASWF010000012.1_ORF_f1_6252_HIGH
#   ./scripts/smoke_gpu.sh

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

BASELINE_IDS="${CASCADE_BASELINE_IDS:-NZ_JAASWF010000012.1_ORF_f1_6252_HIGH}"
export CASCADE_WORKERS="${CASCADE_WORKERS:-1}"
export CASCADE_MAX_GENERATIONS="${CASCADE_MAX_GENERATIONS:-2}"
export CASCADE_VARIANTS_PER_GEN="${CASCADE_VARIANTS_PER_GEN:-2}"
export CASCADE_STAGNATION_LIMIT="${CASCADE_STAGNATION_LIMIT:-2}"

echo "== CASCADE GPU smoke =="
echo "  baselines: $BASELINE_IDS"
echo "  workers=$CASCADE_WORKERS gens=$CASCADE_MAX_GENERATIONS variants/gen=$CASCADE_VARIANTS_PER_GEN"
echo "  EVAL_CMD=${EVAL_CMD:-"(auto-detect)"}"
echo "  PXDESIGN_CMD=${PXDESIGN_CMD:-"(unset)"}"

if [[ -f scripts/cascade_env.sh ]]; then
  # shellcheck disable=SC1091
  source scripts/cascade_env.sh || true
fi

python scripts/evolution_orchestrator.py \
  --baseline-ids "$BASELINE_IDS" \
  --max-generations "$CASCADE_MAX_GENERATIONS" \
  --run-id "smoke_$(date +%Y%m%d_%H%M%S)"

echo "== Smoke finished. Inspect outputs/runs/ and outputs/rl_gym_data/ =="
