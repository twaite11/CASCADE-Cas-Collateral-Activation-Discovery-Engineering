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

## Step 3: Two-Environment Setup (Required for Phase 2)

Phase 2 evolution requires **two separate conda environments**:

| Env | Protenix | Purpose |
|-----|----------|---------|
| **cascade-pxdesign** | v0.5.0+pxd | PXDesign variant generation; Protenix mini (OFF/ON filter) |
| **cascade-protenix** | v1.0.4+ | Protenix base model (high-fidelity ternary scoring); **required** — use latest (pip install protenix) |

PXDesign requires the Protenix 0.5.0+pxd fork. The base model (v1.0.0) is needed for accurate structure scoring. The orchestrator runs in **cascade-protenix** and invokes the other env via env vars.

### 3a. Create cascade-pxdesign (PXDesign + Protenix 0.5.0)

```bash
# Clone and install PXDesign (creates env with Protenix 0.5.0+pxd)
cd /workspace
git clone https://github.com/bytedance/PXDesign.git PXDesign_aa_bias_RL
cd PXDesign_aa_bias_RL
bash -x install.sh --env cascade-pxdesign --pkg_manager conda --cuda-version 12.1

# Install CASCADE Python deps
conda activate cascade-pxdesign
pip install -r /path/to/CASCADE/requirements.txt
```

### 3b. Create cascade-protenix (Protenix latest + base model)

```bash
conda create -n cascade-protenix python=3.11
conda activate cascade-protenix
pip install protenix  # Install latest (v1.0.4+); includes protenix_base_default_v1.0.0
# Download the base model when Protenix prompts on first run (or per Protenix docs)
```

Use `pip install protenix` to get the latest release. Base inference requires Protenix 1.0+ (no fallback to 0.5.0). Optional: set `PROTENIX_BASE_MODEL=protenix_base_20250630_v1.0.0` for the applied model with 2025-06-30 data cutoff (better practical performance).

### 3c. Environment Variables for Evolution

Run the evolution orchestrator in **cascade-protenix** and point PXDesign/mini to **cascade-pxdesign**:

```bash
conda activate cascade-protenix

export PXDESIGN_CMD="/root/miniconda3/envs/cascade-pxdesign/bin/pxdesign"
export PROTENIX_CMD="/root/miniconda3/envs/cascade-pxdesign/bin/protenix"

cd /workspace/CASCADE/scripts
python evolution_orchestrator.py
```

- **PXDESIGN_CMD**: PXDesign runs from cascade-pxdesign (0.5.0+pxd required)
- **PROTENIX_CMD**: Protenix *mini* (OFF/ON filter) runs from cascade-pxdesign
- **Base model**: Uses current env's `protenix` (cascade-protenix, v1.0.4+). Default: `protenix_base_default_v1.0.0`. For better practical performance, set `PROTENIX_BASE_MODEL=protenix_base_20250630_v1.0.0`.

Adjust paths if your conda root differs (e.g. `$CONDA_PREFIX/../cascade-pxdesign/bin/...`).

### 3d. Verify

```bash
# In cascade-pxdesign
conda activate cascade-pxdesign
pxdesign --help
protenix predict -h

# In cascade-protenix
conda activate cascade-protenix
protenix pred -h   # or protenix predict -h
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
source .venv/bin/activate   # or: conda activate pxdesign
cd scripts

# Step 1a: Parse FASTAs, build DB, generate JSONs (CPU only)
python 01_parse_and_annotate.py

# Step 1b: Phase 1 screening - Protenix-mini on each hit (GPU, no MSA for speed)
chmod +x 02_run_screening.sh
./02_run_screening.sh
```

Phase 1 uses `SKIP_MSA=1` by default (~5–10 min per hit). Set `SKIP_MSA=0` to use MSA (slower: ~30–50 min per hit via remote server).

**Option 4 – Two-Phase MSA:** Screen all hits fast (no MSA), then re-run top N with MSA for higher-quality seeds:

```bash
# After 02_run_screening.sh completes, re-run top 5 baselines WITH MSA (~30–50 min each)
chmod +x 02b_rerun_top_with_msa.sh
./02b_rerun_top_with_msa.sh 5   # or: ./02b_rerun_top_with_msa.sh 10
```

Evolution then uses the MSA-enhanced PDBs for those top baselines.

### Phase 2: Evolution Loop (GPU)

```bash
conda activate cascade-protenix
export PXDESIGN_CMD="/root/miniconda3/envs/cascade-pxdesign/bin/pxdesign"
export PROTENIX_CMD="/root/miniconda3/envs/cascade-pxdesign/bin/protenix"

cd scripts
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
# Optional: ./02b_rerun_top_with_msa.sh 5 2>&1 | tee ../logs/phase1b_msa_rerun.log
python evolution_orchestrator.py 2>&1 | tee ../logs/phase2_evolution.log
```

---

## Outputs

| Directory | Contents |
|-----------|----------|
| `outputs/phase1_screening/` | Baseline PDBs from screening |
| `outputs/generation_queue/` | PXDesign variant FASTAs |
| `outputs/optimized_switches/` | Elite Cas13-like sequences & PDBs |
| `outputs/rl_gym_data/` | MPNN bias matrices, `rl_training_dataset.jsonl` for post-training (see RL_TRAINING_FORMAT.md) |

---

## Troubleshooting

| Issue | Fix |
|-------|-----|
| `protenix: command not found` | Ensure cascade-protenix is activated for Phase 2 |
| `Base model not supported` / `requires Protenix 1.0.0+` | Run Phase 2 in cascade-protenix; base inference does not fall back to 0.5.0 |
| OOM on Protenix base / PXDesign | Use A100 80GB; increase `SLEEP_AFTER_*` in evolution_orchestrator |
| `pxdesign: command not found` | Set `PXDESIGN_CMD` to cascade-pxdesign's pxdesign binary |
| `No Phase 1 PDB found` | Run `02_run_screening.sh` first; check `outputs/phase1_screening/` |

---

## Quick Reference

```bash
# One-time setup: create cascade-pxdesign and cascade-protenix (see Step 3)

# Phase 1 (use cascade-pxdesign or venv with protenix)
conda activate cascade-pxdesign  # or: source .venv/bin/activate
cd scripts
python 01_parse_and_annotate.py
./02_run_screening.sh

# Optional: Re-run top 5 with MSA for higher-quality seeds (~2.5–4 h)
./02b_rerun_top_with_msa.sh 5

# Phase 2: Run evolution in cascade-protenix, point to cascade-pxdesign for PXDesign/mini
conda activate cascade-protenix
export PXDESIGN_CMD="/root/miniconda3/envs/cascade-pxdesign/bin/pxdesign"
export PROTENIX_CMD="/root/miniconda3/envs/cascade-pxdesign/bin/protenix"
python evolution_orchestrator.py
```
