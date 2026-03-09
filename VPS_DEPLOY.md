<div align="center">

# 🖥️ CASCADE — VPS Deployment Guide

**Deploy the full CASCADE pipeline on a GPU VPS**

RunPod · Lambda · Vast.ai · Any CUDA 12.x host

</div>

---

## Prerequisites

| Requirement | Minimum | Recommended |
|:------------|:--------|:------------|
| **GPU** | A100 40GB | A100 80GB |
| **OS** | Ubuntu 22.04+ | Ubuntu 22.04 LTS |
| **CUDA** | 12.1 | 12.4 |
| **Python** | 3.11 | 3.11 |
| **Conda** | Miniconda | Miniconda |
| **Storage** | 100 GB | 200 GB |

---

## Why Two Environments?

PXDesign (variant generation) depends on **Protenix 0.5.0+pxd**, while the evaluation pipeline needs **Protenix 1.0.4**. These versions are incompatible in a single environment.

```
┌──────────────────────────────┐     ┌──────────────────────────────┐
│  cascade (conda env)         │     │  pxdesign (conda env)        │
│                              │     │                              │
│  Protenix 1.0.4              │     │  Protenix 0.5.0+pxd          │
│  biopython, numpy            │     │  PXDesign CLI                │
│                              │     │                              │
│  ✅ Structure evaluation     │────▶│  ✅ Variant generation       │
│  ✅ Fitness scoring          │     │                              │
│  ✅ Evolution orchestrator   │     │  Called via PXDESIGN_CMD      │
└──────────────────────────────┘     └──────────────────────────────┘
```

The pipeline runs entirely in `cascade`. When it needs PXDesign, it calls into the `pxdesign` env automatically via `PXDESIGN_CMD`.

---

## Step 1 — Clone

```bash
cd /workspace
git clone https://github.com/twaite11/CASCADE-Cas-Collateral-Activation-Discovery-Engineering.git CASCADE
cd CASCADE
```

---

## Step 2 — Automated Dual-Env Setup

```bash
chmod +x scripts/setup_dual_env.sh
./scripts/setup_dual_env.sh
```

<details>
<summary><b>What this does</b></summary>

1. Creates `cascade` conda env → installs Protenix 1.0.4 + biopython + numpy
2. Clones PXDesign (if not found) → runs its `install.sh` to create `pxdesign` env
3. Creates all output directories
4. Generates `scripts/cascade_env.sh` activation helper

</details>

<details>
<summary><b>Options</b></summary>

```bash
# Only set up one env
./scripts/setup_dual_env.sh --cascade
./scripts/setup_dual_env.sh --pxdesign

# Point to existing PXDesign clone
PXDESIGN_REPO=/path/to/PXDesign ./scripts/setup_dual_env.sh

# Specify CUDA version (default: 12.1)
CUDA_VERSION=12.4 ./scripts/setup_dual_env.sh
```

</details>

---

## Step 3 — Activate (Every Session)

```bash
source scripts/cascade_env.sh
```

This does three things:
1. Activates the `cascade` conda env (Protenix 1.0.4)
2. Exports `PXDESIGN_CMD="conda run --no-banner -n pxdesign pxdesign"`
3. Prints a verification summary

<details>
<summary><b>Manual activation (alternative)</b></summary>

```bash
conda activate cascade
export PXDESIGN_CMD="conda run --no-banner -n pxdesign pxdesign"
```

</details>

**Verify:**

```bash
protenix pred -h              # Protenix 1.0.4 help
$PXDESIGN_CMD --help          # PXDesign (runs in pxdesign env)
```

---

## Step 4 — Prepare Input Data

Place files in `data/mined_hits/`:

| File | Format | Contents |
|:-----|:-------|:--------|
| `deep_hits_*.fasta` | FASTA | Cas13e-like protein sequences (ORFs) |
| `deep_hits_*_metadata.csv` | CSV | `sequence_id`, `repeat_domains`, `sra_accession`, `score` |

---

## Step 5 — Run the Pipeline

### Full Pipeline (Recommended)

```bash
source scripts/cascade_env.sh
mkdir -p logs
./scripts/run_pipeline.sh 2>&1 | tee "logs/cascade_$(date +%Y%m%d_%H%M%S).log"
```

`run_pipeline.sh` auto-detects the `pxdesign` conda env and exports `PXDESIGN_CMD` if not already set.

### Phase by Phase

```bash
source scripts/cascade_env.sh
cd scripts

# Phase 1a: Parse & annotate (CPU only, ~1 min)
python 01_parse_and_annotate.py

# Phase 1b: Protenix-mini screening (GPU, ~5-10 min/hit)
./02_run_screening.sh

# Phase 1c (optional): Re-run top 5 with MSA (~30-50 min/hit)
./02b_rerun_top_with_msa.sh 5

# Phase 2: Evolution loop (GPU, hours)
python evolution_orchestrator.py
```

---

## How `PXDESIGN_CMD` Works

The wrapper `scripts/03_pxdesign_wrapper.py` reads this environment variable:

```python
pxdesign_bin = os.environ.get("PXDESIGN_CMD", "pxdesign")
```

When set to `"conda run --no-banner -n pxdesign pxdesign"`, every PXDesign call executes inside the `pxdesign` conda env (Protenix 0.5.0+pxd), while all other pipeline code stays in `cascade` (Protenix 1.0.4).

**Performance tip:** `conda run` adds ~2s overhead per call. For faster invocation, use the full binary path:

```bash
export PXDESIGN_CMD="/path/to/miniconda3/envs/pxdesign/bin/pxdesign"
```

---

## Outputs

| Directory | Contents |
|:----------|:---------|
| `outputs/phase1_screening/` | Baseline CIF/PDB structures from screening |
| `outputs/generation_queue/` | PXDesign variant FASTAs |
| `outputs/fast_eval/` | Mini-model OFF/ON structural checks |
| `outputs/high_fidelity_scoring/` | Base-model ternary predictions |
| `outputs/rl_gym_data/` | RL bias matrices + `rl_training_dataset.jsonl` |
| `outputs/optimized_switches/` | 🏆 Elite FASTAs + ternary complexes + crRNAs |

---

## Troubleshooting

<details>
<summary><b><code>protenix: command not found</code></b></summary>

Ensure the cascade env is activated:
```bash
source scripts/cascade_env.sh
# or: conda activate cascade
```

</details>

<details>
<summary><b><code>pxdesign: command not found</code></b></summary>

Run the setup script or set `PXDESIGN_CMD` manually:
```bash
./scripts/setup_dual_env.sh --pxdesign
# or:
export PXDESIGN_CMD="conda run --no-banner -n pxdesign pxdesign"
```

</details>

<details>
<summary><b>OOM (Out of Memory) on Protenix / PXDesign</b></summary>

- Use A100 80GB
- Increase `SLEEP_AFTER_PXDESIGN` and `SLEEP_AFTER_PROTENIX` in `evolution_orchestrator.py`
- Run `gc.collect()` and `torch.cuda.empty_cache()` between phases

</details>

<details>
<summary><b>PXDesign uses wrong Protenix version</b></summary>

Verify the pxdesign env has 0.5.x:
```bash
conda run -n pxdesign python -c "import protenix; print(protenix.__version__)"
# Should show 0.5.x
```

</details>

<details>
<summary><b>Protenix eval uses wrong version</b></summary>

Verify the cascade env has 1.0.4+:
```bash
python -c "import protenix; print(protenix.__version__)"
# Should show 1.0.4+
```

</details>

<details>
<summary><b><code>No Phase 1 PDB found</code></b></summary>

Run Phase 1 screening first:
```bash
cd scripts
./02_run_screening.sh
ls ../outputs/phase1_screening/
```

</details>

---

## Quick Reference

```bash
# ── One-time setup ──
./scripts/setup_dual_env.sh

# ── Every session ──
source scripts/cascade_env.sh

# ── Run ──
cd scripts
python 01_parse_and_annotate.py      # Phase 1a: Parse
./02_run_screening.sh                # Phase 1b: Screen
python evolution_orchestrator.py     # Phase 2:  Evolve

# ── Or all at once ──
./scripts/run_pipeline.sh 2>&1 | tee logs/run.log
```

---

## Alternative: Simple Single-Env Setup

If you only need Protenix (no PXDesign), use the simpler venv:

```bash
./scripts/setup_vps.sh
source .venv/bin/activate
```

> ⚠️ PXDesign won't work in this setup. Set `PXDESIGN_CMD` separately if needed.

---

*See [README.md](README.md) for the full project overview and architecture.*
