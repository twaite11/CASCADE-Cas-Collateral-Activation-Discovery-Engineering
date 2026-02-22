#!/bin/bash
# --- CASCADE Full Pipeline Runner ---
# Run all phases with logging. Execute from project root with venv activated.
# Usage: ./scripts/run_pipeline.sh
#        Or: ./scripts/run_pipeline.sh 2>&1 | tee ../logs/cascade_$(date +%Y%m%d_%H%M%S).log
#
# Optional (Option 4 - Two-Phase MSA): Re-run top N baselines with MSA for higher-quality seeds
#   RUN_MSA_RERUN=1 MSA_RERUN_N=5 ./scripts/run_pipeline.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

log_ts() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

cd "$SCRIPT_DIR"

# Check venv
if [ -z "$VIRTUAL_ENV" ] && [ -d "$PROJECT_ROOT/.venv" ]; then
    log_ts "Activating venv..."
    source "$PROJECT_ROOT/.venv/bin/activate"
fi

log_ts "========== CASCADE Pipeline Start =========="

# Phase 1a: Parse & Annotate
log_ts "Phase 1a: Parse and annotate (CPU)..."
python 01_parse_and_annotate.py

# Phase 1b: Screening
log_ts "Phase 1b: Protenix-Mini screening (GPU, no MSA)..."
./02_run_screening.sh

# Phase 1c (optional): Re-run top N with MSA for higher-quality seeds
if [ "${RUN_MSA_RERUN:-0}" = "1" ]; then
    log_ts "Phase 1c: Re-running top ${MSA_RERUN_N:-5} baselines with MSA..."
    chmod +x 02b_rerun_top_with_msa.sh 2>/dev/null || true
    ./02b_rerun_top_with_msa.sh "${MSA_RERUN_N:-5}"
fi

# Phase 2: Evolution loop
log_ts "Phase 2: Evolution loop (GPU)..."
python evolution_orchestrator.py

log_ts "========== CASCADE Pipeline Complete =========="
