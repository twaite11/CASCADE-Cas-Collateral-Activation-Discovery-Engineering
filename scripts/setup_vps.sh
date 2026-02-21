#!/bin/bash
# --- CASCADE VPS Setup Script ---
# Creates a venv, installs dependencies, and optionally configures PXDesign.
# Run from project root: ./scripts/setup_vps.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VENV_DIR="$PROJECT_ROOT/.venv"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] CASCADE VPS Setup"
echo "Project root: $PROJECT_ROOT"

# Ensure we're in project root
cd "$PROJECT_ROOT"

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

echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setup complete."
echo ""
echo "To activate the venv manually:"
echo "  source $VENV_DIR/bin/activate"
echo ""
echo "PXDesign (variant generation) must be installed separately:"
echo "  See VPS_DEPLOY.md for PXDesign installation via conda or Docker."
echo ""
