#!/bin/bash
# --- CASCADE VPS Setup Script ---
# Quick single-env setup using a Python venv.
# For production two-env setup (Protenix 1.0.4 + PXDesign/Protenix 0.5.0), use:
#   ./scripts/setup_dual_env.sh
#
# Run from project root: ./scripts/setup_vps.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VENV_DIR="$PROJECT_ROOT/.venv"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] CASCADE VPS Setup"
echo "Project root: $PROJECT_ROOT"

# Ensure we're in project root
cd "$PROJECT_ROOT"

# Check if conda dual-env setup might be more appropriate
if command -v conda &>/dev/null; then
    echo ""
    echo "NOTE: conda is available. For the recommended dual-env setup"
    echo "      (Protenix 1.0.4 + PXDesign with Protenix 0.5.0), run instead:"
    echo ""
    echo "      ./scripts/setup_dual_env.sh"
    echo ""
    echo "Continuing with simple venv setup..."
    echo ""
fi

# Create venv if it doesn't exist
if [ ! -d "$VENV_DIR" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating virtual environment at $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Virtual environment already exists."
fi

# Activate venv
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Activating venv..."
source "$VENV_DIR/bin/activate"

# Upgrade pip
pip install --upgrade pip -q

# Install requirements
if [ -f "$PROJECT_ROOT/requirements.txt" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Installing dependencies from requirements.txt..."
    pip install -r "$PROJECT_ROOT/requirements.txt"
else
    echo "Warning: requirements.txt not found. Installing core deps manually."
    pip install numpy biopython protenix
fi

# Verify Protenix
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Verifying Protenix..."
if command -v protenix &>/dev/null; then
    echo "  Protenix CLI: OK"
else
    echo "  Protenix CLI: Not in PATH (may need to run from venv)"
fi

# Create output directories
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating output directories..."
mkdir -p "$PROJECT_ROOT/outputs/phase1_screening"
mkdir -p "$PROJECT_ROOT/outputs/generation_queue"
mkdir -p "$PROJECT_ROOT/outputs/fast_eval"
mkdir -p "$PROJECT_ROOT/outputs/high_fidelity_scoring"
mkdir -p "$PROJECT_ROOT/outputs/optimized_switches"
mkdir -p "$PROJECT_ROOT/outputs/rl_gym_data"
mkdir -p "$PROJECT_ROOT/metadata"
mkdir -p "$PROJECT_ROOT/databases"
mkdir -p "$PROJECT_ROOT/logs"

echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setup complete."
echo ""
echo "To activate the venv manually:"
echo "  source $VENV_DIR/bin/activate"
echo ""
echo "PXDesign (variant generation) requires Protenix 0.5.0+pxd in a separate env."
echo "For the recommended dual-env setup:  ./scripts/setup_dual_env.sh"
echo "Or see VPS_DEPLOY.md for manual PXDesign installation."
echo ""
echo "If PXDesign is in a separate conda env, export PXDESIGN_CMD before running:"
echo "  export PXDESIGN_CMD=\"conda run --no-banner -n pxdesign pxdesign\""
echo ""
