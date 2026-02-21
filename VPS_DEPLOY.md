# CASCADE VPS Deployment Guide

Deploy the CASCADE pipeline on a GPU VPS (e.g., RunPod, Lambda, Vast.ai).

---

## Prerequisites

- **GPU**: A100 80GB (recommended) or A100 40GB / L40
- **OS**: Ubuntu 22.04+ with CUDA 12.x
- **Python**: 3.11+
- **Storage**: 100–200 GB free

---

## Step 1: Clone and Enter Project

```bash
cd /workspace   # or your preferred path
git clone <your-repo-url> CASCADE
cd CASCADE
```

---

## Step 2: Create Virtual Environment & Install Dependencies

```bash
# Run the setup script (creates venv, installs deps)
chmod +x scripts/setup_vps.sh
./scripts/setup_vps.sh

# Activate the venv for all subsequent steps
source .venv/bin/activate
```

This installs: `numpy`, `biopython`, `protenix`.

---

## Step 3: Install PXDesign (Required for Phase 2)

PXDesign is used for variant generation. Install via conda (recommended by PXDesign):

```bash
# Clone and install PXDesign
git clone https://github.com/bytedance/PXDesign.git
cd PXDesign
bash -x install.sh --env pxdesign --pkg_manager conda --cuda-version 12.1
conda activate pxdesign

# Install CASCADE Python deps into the same env, or use the venv and ensure pxdesign CLI is available
pip install -r /path/to/CASCADE/requirements.txt
cd /path/to/CASCADE
```

**Alternative**: Use PXDesign's Docker image and call it from your venv (ensure `pxdesign` is on PATH).

Verify:

```bash
protenix pred -h   # Protenix
pxdesign --help    # PXDesign (if installed)
```

---

## Step 4: Prepare Input Data

Place your data in `data/mined_hits/`:

- `deep_hits_*.fasta` — Cas13e-like protein sequences
- `deep_hits_*_metadata.csv` — Metadata with `sequence_id`, `repeat_domains`, `sra_accession`, `score`

Optional: Add MSA databases to `databases/` for better Protenix inputprep (see README).

---

## Step 5: Run the Pipeline

All commands should be run from the project root with the venv activated.

### Phase 1: Parse & Screen (CPU + GPU)

```bash
source .venv/bin/activate
cd scripts

# Step 1a: Parse FASTAs, build DB, generate JSONs (CPU only)
python 01_parse_and_annotate.py

# Step 1b: Phase 1 screening - Protenix-mini on each hit (GPU)
./02_run_screening.sh
```

Phase 1 can take 5–30+ minutes per JSON depending on GPU.

### Phase 2: Evolution Loop (GPU)

```bash
python evolution_orchestrator.py
```

Runs until `MAX_GENERATIONS` or lineage queue is empty. Logs show progress for each step.

---

## Step 6: Run with Logging to File (Recommended)

All steps now emit timestamped logs (e.g. `[2025-02-21 14:30:00] INFO ...`). To capture output:

```bash
mkdir -p logs
source .venv/bin/activate
cd scripts

# Full pipeline with timestamped log
chmod +x run_pipeline.sh
./run_pipeline.sh 2>&1 | tee "../logs/cascade_$(date +%Y%m%d_%H%M%S).log"
```

Or run phases separately:

```bash
python 01_parse_and_annotate.py 2>&1 | tee ../logs/phase1_parse.log
./02_run_screening.sh 2>&1 | tee ../logs/phase1_screen.log
python evolution_orchestrator.py 2>&1 | tee ../logs/phase2_evolution.log
```

---

## Outputs

| Directory | Contents |
|-----------|----------|
| `outputs/phase1_screening/` | Baseline PDBs from screening |
| `outputs/generation_queue/` | PXDesign variant FASTAs |
| `outputs/optimized_switches/` | Elite Cas13-like sequences & PDBs |
| `outputs/rl_gym_data/` | MPNN bias matrices |

---

## Troubleshooting

| Issue | Fix |
|-------|-----|
| `protenix: command not found` | Ensure venv is activated; `pip install protenix` |
| OOM on Protenix base / PXDesign | Use A100 80GB; increase `SLEEP_AFTER_*` in evolution_orchestrator |
| `pxdesign: command not found` | Install PXDesign (Step 3); add to PATH |
| `No Phase 1 PDB found` | Run `02_run_screening.sh` first; check `outputs/phase1_screening/` |

---

## Quick Reference

```bash
# One-time setup
./scripts/setup_vps.sh
source .venv/bin/activate

# Full run
cd scripts
python 01_parse_and_annotate.py
./02_run_screening.sh
python evolution_orchestrator.py
```
