<div align="center">

# 🐳 CASCADE — RunPod Docker Setup

**Build, push, and run CASCADE on RunPod GPU pods**

</div>

---

## Overview

The Docker image packages both environments into a single container:

| Layer | What | Version |
|:------|:-----|:--------|
| **Main** (`/app`) | Protenix + CASCADE pipeline | Protenix 1.0.4 |
| **Isolated** (`/opt/pxdesign_env`) | PXDesign + Protenix | 0.5.0+pxd |

The pipeline runs in the main Python environment. When it needs PXDesign, it calls `/opt/pxdesign_env/bin/pxdesign` via the `PXDESIGN_CMD` environment variable — no conda, no env switching needed inside the container.

---

## Prerequisites

- Docker installed locally
- Docker Hub account (or other container registry)
- RunPod account with GPU pod access

---

## 1. Build the Docker Image

From the project root:

```bash
docker build --platform linux/amd64 -t cascade-runpod:latest .
```

> `--platform linux/amd64` is required when building on ARM (Mac M-series, Windows ARM). RunPod runs x86_64.

Build time: ~15–30 minutes (depends on network for PXDesign deps + JAX CUDA wheels).

---

## 2. Test Locally

```bash
docker run --rm -it --gpus all cascade-runpod:latest /bin/bash
```

Inside the container:

```bash
# Verify main environment (Protenix 1.0.4)
python -c "import protenix; print('Protenix:', protenix.__version__)"
protenix pred -h

# Verify PXDesign environment (Protenix 0.5.0+pxd)
$PXDESIGN_CMD --help

# Verify PXDESIGN_CMD is set
echo $PXDESIGN_CMD
# → /opt/pxdesign_env/bin/pxdesign

# Run tests (no GPU needed)
cd /app && python -m pytest tests/ -v
```

---

## 3. Push to Docker Hub

```bash
docker tag cascade-runpod:latest YOUR_DOCKERHUB_USER/cascade-runpod:latest
docker login
docker push YOUR_DOCKERHUB_USER/cascade-runpod:latest
```

---

## 4. Create a RunPod Template

1. Go to [RunPod Templates](https://console.runpod.io/user/templates) → **New Template**
2. Configure:

| Setting | Value |
|:--------|:------|
| **Name** | `cascade-gpu` |
| **Container Image** | `YOUR_DOCKERHUB_USER/cascade-runpod:latest` |
| **Container Disk** | 50+ GB (for weights + outputs) |
| **TCP Ports** | 22 (SSH) |
| **HTTP Ports** | 8888 (JupyterLab) |

3. Save the template.

---

## 5. First-Time Setup: Download Weights

On the first pod launch, download model weights:

```bash
cd /opt/PXDesign
bash download_tool_weights.sh
```

This downloads:
- AlphaFold2 parameters (~4 GB)
- ProteinMPNN weights
- CCD cache (components.cif, etc.)

> **Tip**: Use a RunPod **Network Volume** to persist weights across pods (see [section 8](#8-network-volume-for-weights-recommended)).

---

## 6. Run the Pipeline

```bash
cd /app/scripts

# Phase 1: Parse & Screen
python 01_parse_and_annotate.py
./02_run_screening.sh

# Phase 2: Evolution
python evolution_orchestrator.py
```

Or all at once with logging:

```bash
cd /app
mkdir -p logs
./scripts/run_pipeline.sh 2>&1 | tee "logs/cascade_$(date +%Y%m%d_%H%M%S).log"
```

---

## 7. Environment Variables

These are pre-configured in the Dockerfile. Override via RunPod template env vars or at runtime.

### Pre-Set (usually no changes needed)

| Variable | Default | Description |
|:---------|:--------|:------------|
| `PXDESIGN_CMD` | `/opt/pxdesign_env/bin/pxdesign` | PXDesign binary path (isolated env) |
| `PROTENIX_DATA_ROOT_DIR` | `/opt/PXDesign/release_data/ccd_cache` | CCD cache path |
| `TOOL_WEIGHTS_ROOT` | `/opt/PXDesign/tool_weights` | AF2, MPNN weight path |
| `CUTLASS_PATH` | `/opt/cutlass` | CUTLASS for DeepSpeed Evo attention |

### Optional Overrides

| Variable | Default | Description |
|:---------|:--------|:------------|
| `PROTENIX_BASE_MODEL` | `protenix_base_default_v1.0.0` | Override Protenix base model |
| `CASCADE_RUN_ID` | timestamp | Run identifier for output folders |

---

## 8. Network Volume for Weights (Recommended)

Avoid re-downloading ~4 GB of weights on every pod:

1. Create a **Network Volume** in RunPod (100 GB)
2. On first pod start with the volume attached:

```bash
# Volume is mounted at /runpod-volume
mkdir -p /runpod-volume/cascade_weights
cd /opt/PXDesign
bash download_tool_weights.sh

# Copy weights to persistent volume
cp -r /opt/PXDesign/tool_weights/* /runpod-volume/cascade_weights/
cp -r /opt/PXDesign/release_data/ccd_cache/* /runpod-volume/cascade_weights/
```

3. Set env vars to point at the volume:

```bash
export TOOL_WEIGHTS_ROOT=/runpod-volume/cascade_weights
export PROTENIX_DATA_ROOT_DIR=/runpod-volume/cascade_weights
```

4. On future pods, just attach the same volume — no download needed.

---

## 9. Quick Reference

```bash
# Build
docker build --platform linux/amd64 -t cascade-runpod:latest .

# Push
docker push YOUR_USER/cascade-runpod:latest

# First pod: download weights
cd /opt/PXDesign && bash download_tool_weights.sh

# Run CASCADE
cd /app/scripts && python evolution_orchestrator.py
```

---

## Troubleshooting

<details>
<summary><b>Out of disk space</b></summary>

Increase container disk to 50+ GB, or use a network volume for weights and outputs.

</details>

<details>
<summary><b>CUDA / GPU errors</b></summary>

Use A100 40GB or 80GB pods. Avoid older GPUs (< Ampere). Check:
```bash
python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
```

</details>

<details>
<summary><b>PXDesign not found</b></summary>

Verify the isolated env:
```bash
$PXDESIGN_CMD --help
# or directly:
/opt/pxdesign_env/bin/pxdesign --help
```

</details>

<details>
<summary><b>Protenix version mismatch</b></summary>

```bash
# Main env should be 1.0.4+
python -c "import protenix; print(protenix.__version__)"

# PXDesign env should be 0.5.x
/opt/pxdesign_env/bin/python -c "import protenix; print(protenix.__version__)"
```

</details>

<details>
<summary><b>Weights not downloaded</b></summary>

```bash
cd /opt/PXDesign && bash download_tool_weights.sh
```

Or set `TOOL_WEIGHTS_ROOT` and `PROTENIX_DATA_ROOT_DIR` to a pre-downloaded volume.

</details>

---

*See [VPS_DEPLOY.md](VPS_DEPLOY.md) for bare-metal conda setup. See [README.md](README.md) for project overview.*
