# CASCADE RunPod Docker Setup Guide

This guide walks you through building the CASCADE Docker image for RunPod and configuring environment variables.

---

## Prerequisites

- Docker installed
- Docker Hub account (or another container registry)
- RunPod account
- **Before building**: Ensure the PXDesign submodule is populated:
  ```bash
  git submodule update --init --recursive
  ```

---

## 1. Build the Docker Image

From the project root:

```bash
docker build --platform linux/amd64 -t cascade-runpod:latest .
```

> **Note**: `--platform linux/amd64` is required if building on Mac (ARM) or Windows ARM; RunPod runs on x86_64.

Build time is typically 15–30 minutes depending on network speed.

---

## 2. Test Locally

```bash
docker run --rm -it --gpus all cascade-runpod:latest /bin/bash
```

Inside the container:

```bash
# Verify installations
python -c "import torch; print('PyTorch:', torch.__version__, '| CUDA:', torch.cuda.is_available())"
python -c "import protenix, pxdesign; print('Protenix & PXDesign OK')"
```

---

## 3. Push to Docker Hub

```bash
# Tag with your Docker Hub username
docker tag cascade-runpod:latest YOUR_DOCKERHUB_USERNAME/cascade-runpod:latest

# Log in
docker login

# Push
docker push YOUR_DOCKERHUB_USERNAME/cascade-runpod:latest
```

---

## 4. Create a RunPod Template

1. Go to [RunPod Templates](https://console.runpod.io/user/templates) → **New Template**
2. Configure:
   - **Name**: `cascade-pxdesign` (or similar)
   - **Container Image**: `YOUR_DOCKERHUB_USERNAME/cascade-runpod:latest`
   - **Container Disk**: **至少 50 GB** (for weights + outputs)
   - **TCP Ports**: SSH = 22
   - **HTTP Ports**: JupyterLab = 8888
3. Save the template.

---

## 5. First-Time Pod Setup: Download Weights

On the **first launch** of a Pod, download model weights inside the container:

```bash
cd /app/PXDesign_aa_bias_RL
bash download_tool_weights.sh
```

This downloads:

- AlphaFold2 parameters (~4 GB)
- ProteinMPNN weights
- CCD cache (components.cif, etc.)

**Tip**: Create a RunPod network volume, run this once, then attach the volume to future Pods to reuse the weights.

---

## 6. Environment Variables

Set these in RunPod **Template → Environment Variables** or at runtime.

### Required for CASCADE

| Variable | Description | Default |
|----------|-------------|---------|
| `PROTENIX_DATA_ROOT_DIR` | Path to CCD cache (components.cif, etc.) | `/app/PXDesign_aa_bias_RL/release_data/ccd_cache` |
| `TOOL_WEIGHTS_ROOT` | Path to AF2, MPNN, etc. | `/app/PXDesign_aa_bias_RL/tool_weights` |

### Optional – Evolution Orchestrator

| Variable | Description | Default |
|----------|-------------|---------|
| `CASCADE_RUN_ID` | Run identifier for output folders | `YYYYMMDD_HHMMSS` |
| `CASCADE_NUM_WORKERS` | Number of parallel lineage workers | Number of visible GPUs |

### Optional – Protenix CLI

| Variable | Description | Default |
|----------|-------------|---------|
| `PROTENIX_CMD` | Protenix executable path | `protenix` |
| `PROTENIX_BASE_MODEL` | Protenix base model name | `protenix_base_default_v1.0.0` |

### Optional – PXDesign

| Variable | Description | Default |
|----------|-------------|---------|
| `PXDESIGN_CMD` | PXDesign executable path | `pxdesign` |
| `CUTLASS_PATH` | CUTLASS directory (DeepSpeed Evo) | `/opt/cutlass` |

---

## 7. RunPod Template Environment Variables (Summary)

For a typical setup, you can leave defaults and only override if needed:

```
# Optional – only if using custom paths
PROTENIX_DATA_ROOT_DIR=/app/PXDesign_aa_bias_RL/release_data/ccd_cache
TOOL_WEIGHTS_ROOT=/app/PXDesign_aa_bias_RL/tool_weights

# Optional – multi-GPU
CASCADE_NUM_WORKERS=4

# Optional – Run ID
CASCADE_RUN_ID=my_experiment_001
```

---

## 8. Using a Network Volume for Weights (Recommended)

1. Create a **Network Volume** in RunPod (e.g. 100 GB).
2. On first Pod start with the volume attached:
   ```bash
   # Volume is mounted at /runpod-volume
   cd /app/PXDesign_aa_bias_RL
   mkdir -p /runpod-volume/cascade_weights
   bash download_tool_weights.sh /runpod-volume/cascade_weights
   ```
   This puts AF2, MPNN, and CCD cache under `/runpod-volume/cascade_weights/`.
3. Set environment variables to point to the volume:
   ```
   TOOL_WEIGHTS_ROOT=/runpod-volume/cascade_weights
   PROTENIX_DATA_ROOT_DIR=/runpod-volume/cascade_weights
   ```
4. Attach the same volume to future Pods to skip re-downloading.

---

## 9. Running the Evolution Pipeline

Once weights are downloaded:

```bash
cd /app
python scripts/evolution_orchestrator.py
```

Or with custom run ID and worker count:

```bash
CASCADE_RUN_ID=exp_001 CASCADE_NUM_WORKERS=4 python scripts/evolution_orchestrator.py
```

---

## 10. Quick Reference

| Step | Action |
|------|--------|
| Build | `docker build --platform linux/amd64 -t cascade-runpod:latest .` |
| Push | `docker push YOUR_USER/cascade-runpod:latest` |
| First run | `cd /app/PXDesign_aa_bias_RL && bash download_tool_weights.sh` |
| Run CASCADE | `python scripts/evolution_orchestrator.py` |

---

## Troubleshooting

- **Out of disk space**: Increase container disk to 50+ GB or use a network volume.
- **CUDA errors**: Choose a Pod with an A100 or similar GPU; avoid very old GPUs.
- **Protenix not found**: Ensure `download_tool_weights.sh` completed and paths are correct.
- **Submodule empty on build**: Run `git submodule update --init --recursive` before `docker build`.
