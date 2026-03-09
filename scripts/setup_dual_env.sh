#!/bin/bash
# =============================================================================
# CASCADE Dual-Environment Setup Script
# =============================================================================
# Creates TWO isolated conda environments:
#
#   1. "cascade"  — Main pipeline env (Protenix 1.0.4 for structure evaluation)
#   2. "pxdesign" — PXDesign env (Protenix 0.5.0+pxd for variant generation)
#
# The evolution orchestrator runs inside "cascade" and calls PXDesign from
# the "pxdesign" env via:  PXDESIGN_CMD="conda run --no-banner -n pxdesign pxdesign"
#
# Usage:
#   chmod +x scripts/setup_dual_env.sh
#   ./scripts/setup_dual_env.sh              # full setup (both envs)
#   ./scripts/setup_dual_env.sh --cascade    # cascade env only
#   ./scripts/setup_dual_env.sh --pxdesign   # pxdesign env only
#
# After setup, source the activation helper:
#   source scripts/cascade_env.sh
# =============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

log_ts() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
SETUP_CASCADE=true
SETUP_PXDESIGN=true
PXDESIGN_REPO="${PXDESIGN_REPO:-}"  # Optional: path to existing PXDesign clone

if [ "$1" = "--cascade" ]; then
    SETUP_PXDESIGN=false
elif [ "$1" = "--pxdesign" ]; then
    SETUP_CASCADE=false
fi

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
log_ts "${CYAN}CASCADE Dual-Environment Setup${NC}"
log_ts "Project root: $PROJECT_ROOT"

if ! command -v conda &>/dev/null; then
    log_ts "${RED}ERROR: conda not found. Install Miniconda or Anaconda first:${NC}"
    echo "  https://docs.anaconda.com/miniconda/"
    exit 1
fi

# Ensure conda shell hooks are available
eval "$(conda shell.bash hook 2>/dev/null)" || true

CUDA_VERSION="${CUDA_VERSION:-12.1}"
PYTHON_VERSION="${PYTHON_VERSION:-3.11}"

log_ts "CUDA version target: ${CUDA_VERSION}"
log_ts "Python version: ${PYTHON_VERSION}"

# ---------------------------------------------------------------------------
# 1. CASCADE environment (Protenix 1.0.4 — structure evaluation)
# ---------------------------------------------------------------------------
if [ "$SETUP_CASCADE" = true ]; then
    log_ts ""
    log_ts "${GREEN}=== Setting up 'cascade' environment (Protenix 1.0.4) ===${NC}"

    if conda env list | grep -q "^cascade "; then
        log_ts "${YELLOW}conda env 'cascade' already exists. Updating...${NC}"
        conda activate cascade
    else
        log_ts "Creating conda env 'cascade' (Python ${PYTHON_VERSION})..."
        conda create -n cascade python="${PYTHON_VERSION}" -y -q
        conda activate cascade
    fi

    log_ts "Installing Protenix 1.0.4 and CASCADE dependencies..."
    pip install --upgrade pip -q
    pip install "protenix>=1.0.4,<2.0" -q
    pip install -r "$PROJECT_ROOT/requirements.txt" -q

    # Verify
    PROTENIX_VER=$(python -c "import protenix; print(protenix.__version__)" 2>/dev/null || echo "unknown")
    log_ts "Protenix version in cascade env: ${GREEN}${PROTENIX_VER}${NC}"

    if command -v protenix &>/dev/null; then
        log_ts "protenix CLI: ${GREEN}OK${NC}"
    else
        log_ts "${YELLOW}protenix CLI not on PATH (may work via 'python -m protenix')${NC}"
    fi

    conda deactivate
    log_ts "${GREEN}'cascade' environment ready.${NC}"
fi

# ---------------------------------------------------------------------------
# 2. PXDESIGN environment (Protenix 0.5.0+pxd — variant generation)
# ---------------------------------------------------------------------------
if [ "$SETUP_PXDESIGN" = true ]; then
    log_ts ""
    log_ts "${GREEN}=== Setting up 'pxdesign' environment (Protenix 0.5.0+pxd) ===${NC}"

    # Check if pxdesign env already exists
    if conda env list | grep -q "^pxdesign "; then
        log_ts "${YELLOW}conda env 'pxdesign' already exists. Skipping creation.${NC}"
        log_ts "To recreate: conda env remove -n pxdesign && rerun this script"
    else
        # PXDesign has its own install.sh that creates the conda env
        if [ -n "$PXDESIGN_REPO" ] && [ -d "$PXDESIGN_REPO" ]; then
            PXDESIGN_DIR="$PXDESIGN_REPO"
            log_ts "Using existing PXDesign clone: $PXDESIGN_DIR"
        elif [ -d "$PROJECT_ROOT/../PXDesign" ]; then
            PXDESIGN_DIR="$(cd "$PROJECT_ROOT/../PXDesign" && pwd)"
            log_ts "Found PXDesign at: $PXDESIGN_DIR"
        else
            log_ts "Cloning PXDesign repository..."
            PXDESIGN_DIR="$PROJECT_ROOT/../PXDesign"
            git clone https://github.com/bytedance/PXDesign.git "$PXDESIGN_DIR"
        fi

        log_ts "Running PXDesign install.sh (creates 'pxdesign' conda env with Protenix 0.5.0+pxd)..."
        log_ts "${YELLOW}This may take 10-20 minutes...${NC}"
        cd "$PXDESIGN_DIR"
        bash -x install.sh --env pxdesign --pkg_manager conda --cuda-version "$CUDA_VERSION"
        cd "$PROJECT_ROOT"
    fi

    # Verify pxdesign CLI is available in the env
    log_ts "Verifying pxdesign CLI..."
    if conda run --no-banner -n pxdesign pxdesign --help &>/dev/null; then
        log_ts "pxdesign CLI: ${GREEN}OK${NC}"
    else
        log_ts "${YELLOW}pxdesign CLI not responding via conda run. You may need to check the install.${NC}"
    fi

    # Verify Protenix version in pxdesign env
    PXD_PROTENIX_VER=$(conda run --no-banner -n pxdesign python -c \
        "import protenix; print(protenix.__version__)" 2>/dev/null || echo "unknown")
    log_ts "Protenix version in pxdesign env: ${GREEN}${PXD_PROTENIX_VER}${NC}"

    log_ts "${GREEN}'pxdesign' environment ready.${NC}"
fi

# ---------------------------------------------------------------------------
# 3. Create output directories
# ---------------------------------------------------------------------------
log_ts ""
log_ts "Creating output directories..."
mkdir -p "$PROJECT_ROOT/outputs/phase1_screening"
mkdir -p "$PROJECT_ROOT/outputs/generation_queue"
mkdir -p "$PROJECT_ROOT/outputs/fast_eval"
mkdir -p "$PROJECT_ROOT/outputs/high_fidelity_scoring"
mkdir -p "$PROJECT_ROOT/outputs/optimized_switches"
mkdir -p "$PROJECT_ROOT/outputs/rl_gym_data"
mkdir -p "$PROJECT_ROOT/metadata"
mkdir -p "$PROJECT_ROOT/databases"
mkdir -p "$PROJECT_ROOT/logs"

# ---------------------------------------------------------------------------
# 4. Generate the activation helper script
# ---------------------------------------------------------------------------
ACTIVATE_SCRIPT="$SCRIPT_DIR/cascade_env.sh"
log_ts "Writing activation helper: $ACTIVATE_SCRIPT"

cat > "$ACTIVATE_SCRIPT" << 'ENVEOF'
#!/bin/bash
# =============================================================================
# CASCADE Environment Activation Helper
# =============================================================================
# Source this file to activate the cascade env and configure PXDESIGN_CMD:
#
#   source scripts/cascade_env.sh
#
# This sets up:
#   - Activates the 'cascade' conda env (Protenix 1.0.4)
#   - Exports PXDESIGN_CMD so the pipeline calls PXDesign from the 'pxdesign' env
#   - Exports PROJECT_ROOT for convenience
# =============================================================================

_CASCADE_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_ROOT="$(cd "$_CASCADE_SCRIPT_DIR/.." && pwd)"

# Activate cascade conda env
eval "$(conda shell.bash hook 2>/dev/null)" || true
conda activate cascade 2>/dev/null

if [ "$CONDA_DEFAULT_ENV" != "cascade" ]; then
    echo "[WARNING] Failed to activate 'cascade' conda env. Is it installed?"
    echo "  Run: ./scripts/setup_dual_env.sh"
    return 1 2>/dev/null || exit 1
fi

# Configure PXDESIGN_CMD to call pxdesign from its isolated conda env
# This is read by scripts/03_pxdesign_wrapper.py (line: os.environ.get("PXDESIGN_CMD"))
export PXDESIGN_CMD="conda run --no-banner -n pxdesign pxdesign"

# Verify both tools are reachable
echo "============================================="
echo " CASCADE Environment Active"
echo "============================================="
echo " conda env:    $CONDA_DEFAULT_ENV"

_PROT_VER=$(python -c "import protenix; print(protenix.__version__)" 2>/dev/null || echo "?")
echo " protenix:     v${_PROT_VER} (evaluation)"

_PXD_VER=$(conda run --no-banner -n pxdesign pxdesign --version 2>/dev/null || echo "?")
echo " pxdesign:     ${_PXD_VER} (generation, via pxdesign env)"
echo " PXDESIGN_CMD: $PXDESIGN_CMD"
echo " PROJECT_ROOT: $PROJECT_ROOT"
echo "============================================="
echo ""
echo "Ready. Run the pipeline with:"
echo "  cd scripts && python evolution_orchestrator.py"
echo "  # or: ./scripts/run_pipeline.sh"

unset _CASCADE_SCRIPT_DIR _PROT_VER _PXD_VER
ENVEOF

chmod +x "$ACTIVATE_SCRIPT"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log_ts ""
log_ts "${GREEN}=============================================${NC}"
log_ts "${GREEN} Setup Complete!${NC}"
log_ts "${GREEN}=============================================${NC}"
echo ""
echo "  Two conda environments are now configured:"
echo ""
echo "    cascade   — Protenix 1.0.4 (structure eval, fitness scoring)"
echo "    pxdesign  — Protenix 0.5.0+pxd (PXDesign variant generation)"
echo ""
echo "  The pipeline runs in 'cascade' and calls PXDesign cross-env."
echo ""
echo "  To activate for every session:"
echo "    source scripts/cascade_env.sh"
echo ""
echo "  Or manually:"
echo "    conda activate cascade"
echo "    export PXDESIGN_CMD=\"conda run --no-banner -n pxdesign pxdesign\""
echo ""

