# =============================================================================
# CASCADE RunPod Docker Image
# =============================================================================
# Two-stage install:
#   1. Protenix 1.0.4  — structure evaluation (main pipeline)
#   2. PXDesign         — variant generation (stock, from upstream repo)
#
# PXDesign needs Protenix 0.5.0+pxd; we install it into a separate virtualenv
# and call it via PXDESIGN_CMD from the main pipeline.
#
# Build:  docker build --platform linux/amd64 -t cascade-runpod:latest .
# Run:    docker run --gpus all -it cascade-runpod:latest
# =============================================================================

# RunPod PyTorch base (CUDA 12.4, Ubuntu 22.04) — includes Jupyter, SSH
FROM runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04

ENV PYTHONUNBUFFERED=1
WORKDIR /app

# ---------------------------------------------------------------------------
# System dependencies
# ---------------------------------------------------------------------------
RUN apt-get update --yes && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
        git \
        curl \
        wget \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 1. CASCADE main environment (Protenix 1.0.4)
# ---------------------------------------------------------------------------
COPY requirements.txt /app/
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir "protenix>=1.0.4,<2.0"

# Verify Protenix 1.0.4
RUN python -c "import protenix; print(f'Protenix {protenix.__version__} installed')"

# ---------------------------------------------------------------------------
# 2. PXDesign environment (separate virtualenv with Protenix 0.5.0+pxd)
# ---------------------------------------------------------------------------
# Create isolated virtualenv for PXDesign so it doesn't conflict with main env
RUN python -m venv /opt/pxdesign_env

# Install Protenix 0.5.0+pxd (PXDesign-compatible) into the isolated env
RUN /opt/pxdesign_env/bin/pip install --no-cache-dir --upgrade pip && \
    /opt/pxdesign_env/bin/pip install --no-cache-dir \
        "git+https://github.com/bytedance/Protenix.git@v0.5.0+pxd"

# PXDesign dependencies
RUN /opt/pxdesign_env/bin/pip install --no-cache-dir \
    einops natsort dm-tree posix_ipc \
    "transformers==4.51.3" "dm-haiku==0.0.13" "optax==0.2.5"

# ColabDesign (no deps)
RUN /opt/pxdesign_env/bin/pip install --no-cache-dir \
    git+https://github.com/sokrypton/ColabDesign.git --no-deps

# JAX with CUDA
RUN /opt/pxdesign_env/bin/pip install --no-cache-dir "jax[cuda]==0.4.29" \
    -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# PXDesignBench
RUN /opt/pxdesign_env/bin/pip install --no-cache-dir "numpy==1.26.3" && \
    /opt/pxdesign_env/bin/pip install --no-cache-dir \
        git+https://github.com/bytedance/PXDesignBench.git@v0.1.2 --no-deps

# Clone and install stock PXDesign (upstream, no custom fork needed)
RUN git clone --depth 1 https://github.com/bytedance/PXDesign.git /opt/PXDesign && \
    /opt/pxdesign_env/bin/pip install --no-cache-dir -e /opt/PXDesign

# Verify PXDesign
RUN /opt/pxdesign_env/bin/pxdesign --help > /dev/null 2>&1 || \
    echo "WARNING: pxdesign CLI not available (may need weights first)"

# CUTLASS (for DeepSpeed Evo attention)
RUN git clone -b v3.5.1 --depth 1 https://github.com/NVIDIA/cutlass.git /opt/cutlass
ENV CUTLASS_PATH=/opt/cutlass

# ---------------------------------------------------------------------------
# 3. Copy CASCADE project
# ---------------------------------------------------------------------------
COPY . /app/

# ---------------------------------------------------------------------------
# 4. Environment configuration
# ---------------------------------------------------------------------------
# PXDESIGN_CMD tells 03_pxdesign_wrapper.py to use the isolated PXDesign env
ENV PXDESIGN_CMD="/opt/pxdesign_env/bin/pxdesign"

# Default data/weight paths (override via RunPod env vars or volume mounts)
ENV PROTENIX_DATA_ROOT_DIR=/opt/PXDesign/release_data/ccd_cache
ENV TOOL_WEIGHTS_ROOT=/opt/PXDesign/tool_weights

# Create directories for weights (download on first use or mount pre-downloaded)
RUN mkdir -p /opt/PXDesign/release_data/ccd_cache \
             /opt/PXDesign/tool_weights \
             /app/outputs \
             /app/logs \
             /app/metadata \
             /app/databases

# Ensure scripts are executable
RUN chmod +x /app/scripts/*.sh 2>/dev/null || true

# Default: keep RunPod's Jupyter + SSH (no CMD override)
# For headless mode: CMD ["python", "/app/scripts/evolution_orchestrator.py"]
