#!/bin/bash
# --- CASCADE Full Pipeline Runner ---
# Run all phases with logging. Execute from project root with env activated.
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

# ---------------------------------------------------------------------------
# Environment activation — supports both venv and conda dual-env setups
# ---------------------------------------------------------------------------
if [ -n "$CONDA_DEFAULT_ENV" ] && [ "$CONDA_DEFAULT_ENV" = "cascade" ]; then
    log_ts "Using active conda env: cascade"
elif [ -z "$VIRTUAL_ENV" ] && [ -d "$PROJECT_ROOT/.venv" ]; then
    log_ts "Activating venv..."
    source "$PROJECT_ROOT/.venv/bin/activate"
elif [ -z "$VIRTUAL_ENV" ] && [ -z "$CONDA_DEFAULT_ENV" ]; then
    # Try cascade conda env
    eval "$(conda shell.bash hook 2>/dev/null)" || true
    if conda env list 2>/dev/null | grep -q "^cascade "; then
        log_ts "Activating 'cascade' conda env..."
        conda activate cascade
    fi
fi

# ---------------------------------------------------------------------------
# Auto-detect PXDESIGN_CMD if not already set
# ---------------------------------------------------------------------------
# PXDesign needs Protenix 0.5.0+pxd (in the 'pxdesign' conda env), while
# the main pipeline uses Protenix 1.0.4 (in 'cascade' or .venv).
# If PXDESIGN_CMD is not set, check if a 'pxdesign' conda env exists and
# configure cross-env invocation automatically.
if [ -z "$PXDESIGN_CMD" ]; then
    if conda env list 2>/dev/null | grep -q "^pxdesign "; then
        export PXDESIGN_CMD="conda run --no-banner -n pxdesign pxdesign"
        log_ts "Auto-detected pxdesign conda env -> PXDESIGN_CMD=$PXDESIGN_CMD"
    elif command -v pxdesign &>/dev/null; then
        log_ts "pxdesign found on PATH (same env)."
    else
        log_ts "WARNING: pxdesign not found. Phase 2 variant generation will fail."
        log_ts "  Run ./scripts/setup_dual_env.sh to install PXDesign, or set PXDESIGN_CMD manually."
    fi
else
    log_ts "Using PXDESIGN_CMD=$PXDESIGN_CMD"
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
