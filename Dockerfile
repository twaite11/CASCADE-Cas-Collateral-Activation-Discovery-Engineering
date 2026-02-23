# CASCADE RunPod Template
# Build: docker build --platform linux/amd64 -t cascade-runpod:latest .
# Run: docker run --gpus all -it cascade-runpod:latest

# RunPod PyTorch base (CUDA 12.4, Ubuntu 22.04) - includes Jupyter, SSH
FROM runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04

ENV PYTHONUNBUFFERED=1
WORKDIR /app

# System dependencies
RUN apt-get update --yes && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
        git \
        curl \
        wget \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

# --- CASCADE core dependencies ---
COPY requirements.txt /app/
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# --- Copy project (required before PXDesign install) ---
# Note: Ensure PXDesign_aa_bias_RL is populated: run `git submodule update --init --recursive` before building
COPY . /app/

# --- PXDesign dependencies (must match PXDesign_aa_bias_RL/install.sh) ---
# Protenix v0.5.0+pxd (PXDesign-compatible fork)
RUN pip install --no-cache-dir "git+https://github.com/bytedance/Protenix.git@v0.5.0+pxd"

# PXDesignBench base deps
RUN pip install --no-cache-dir \
    einops natsort dm-tree posix_ipc \
    "transformers==4.51.3" "dm-haiku==0.0.13" "optax==0.2.5"

# ColabDesign (no deps)
RUN pip install --no-cache-dir git+https://github.com/sokrypton/ColabDesign.git --no-deps

# JAX with CUDA
RUN pip install --no-cache-dir "jax[cuda]==0.4.29" \
    -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# PXDesignBench
RUN pip install --no-cache-dir "numpy==1.26.3" && \
    pip install --no-cache-dir git+https://github.com/bytedance/PXDesignBench.git@v0.1.2 --no-deps

# PXDesign package (from submodule)
RUN pip install --no-cache-dir -e /app/PXDesign_aa_bias_RL

# CUTLASS (for DeepSpeed Evo attention)
RUN git clone -b v3.5.1 --depth 1 https://github.com/NVIDIA/cutlass.git /opt/cutlass
ENV CUTLASS_PATH=/opt/cutlass

# Default data dirs (user can override via env or mount volumes)
ENV PROTENIX_DATA_ROOT_DIR=/app/PXDesign_aa_bias_RL/release_data/ccd_cache
ENV TOOL_WEIGHTS_ROOT=/app/PXDesign_aa_bias_RL/tool_weights

# Create dirs for weights (run download_tool_weights.sh on first use or mount pre-downloaded)
RUN mkdir -p /app/PXDesign_aa_bias_RL/release_data/ccd_cache \
             /app/PXDesign_aa_bias_RL/tool_weights \
             /app/outputs

# Ensure scripts are executable
RUN chmod +x /app/PXDesign_aa_bias_RL/download_tool_weights.sh 2>/dev/null || true

# Default: keep RunPod's Jupyter + SSH (no CMD override)
# For app-only mode, use: CMD ["python", "/app/scripts/evolution_orchestrator.py"]
