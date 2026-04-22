#!/usr/bin/env bash
# =============================================================================
# Orchestrator container entrypoint
# =============================================================================
# Activates the cascade conda env, exports PXDESIGN_CMD so pipeline code can
# call into the pxdesign env, and execs evolution_orchestrator.py with any
# CLI args passed through (or env vars honored by the orchestrator argparse).
#
# Config (all optional, read from env; CLI flags take precedence):
#   CASCADE_RUN_ID            default: "local"
#   CASCADE_BASELINE_IDS      comma-separated
#   CASCADE_CRRNA_LOOKUP_IDS  comma-separated
#   CASCADE_MAX_GENERATIONS   default: 12
#   CASCADE_VARIANTS_PER_GEN  default: 5
#   CASCADE_WORKERS           default: 3
# =============================================================================
set -euo pipefail

# shellcheck disable=SC1091
source /opt/conda/etc/profile.d/conda.sh
conda activate cascade

export PXDESIGN_CMD="${PXDESIGN_CMD:-/opt/conda/envs/pxdesign/bin/pxdesign}"
export CASCADE_ROOT="${CASCADE_ROOT:-/workspace/CASCADE}"

cd "${CASCADE_ROOT}/scripts"

RUN_ID="${CASCADE_RUN_ID:-local}"
LOG_DIR="${CASCADE_ROOT}/logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/run-${RUN_ID}.log"

echo "[entrypoint] cascade env: $(python --version)"
echo "[entrypoint] PXDESIGN_CMD=${PXDESIGN_CMD}"
echo "[entrypoint] run_id=${RUN_ID} logs=${LOG_FILE}"
echo "[entrypoint] exec: python evolution_orchestrator.py $*"

# Tee stdout/stderr to a log file the controller can rsync back, while keeping
# live stream on the container's stdout for the dashboard WebSocket tailer.
exec python -u evolution_orchestrator.py "$@" 2>&1 | tee "${LOG_FILE}"
