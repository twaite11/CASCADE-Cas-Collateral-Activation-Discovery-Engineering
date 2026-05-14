# =============================================================================
# CASCADE Evolution Orchestrator — GPU container image
# =============================================================================
# Builds a self-contained image that runs the directed-evolution orchestrator
# on a Vast.ai A100 VPS (or any CUDA 12.1+ host).
#
# Layout inside the image:
#   /opt/conda               Miniconda root
#   /opt/conda/envs/cascade  Protenix 1.0.4 + CASCADE deps
#   /opt/conda/envs/pxdesign PXDesign + Protenix 0.5.0+pxd
#   /workspace/CASCADE       This repository
#   /workspace/CASCADE/outputs/runs/<RUN_ID>/... per-run artifact tree
#
# Run via Vast.ai by passing --onstart-cmd "..." that invokes
#   scripts/evolution_orchestrator.py with the run-specific flags. A thin
#   orchestrator/entrypoint.sh is available for manual docker run usage too.
# =============================================================================
FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    CONDA_DIR=/opt/conda \
    PATH=/opt/conda/bin:/opt/conda/condabin:$PATH \
    CASCADE_ROOT=/workspace/CASCADE \
    PXDESIGN_CMD=/opt/conda/envs/pxdesign/bin/pxdesign

# -----------------------------------------------------------------------------
# System packages: git+build tools for Python wheels, rsync+openssh for the
# controller's artifact sync, ca-certificates for HTTPS clones.
# -----------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
        wget \
        curl \
        bzip2 \
        rsync \
        openssh-client \
        tini \
        build-essential \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# Miniconda
# -----------------------------------------------------------------------------
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p ${CONDA_DIR} \
    && rm /tmp/miniconda.sh \
    && ${CONDA_DIR}/bin/conda init bash \
    && ${CONDA_DIR}/bin/conda config --set always_yes yes --set changeps1 no \
    && ${CONDA_DIR}/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
    && ${CONDA_DIR}/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r \
    && ${CONDA_DIR}/bin/conda config --add channels conda-forge \
    && ${CONDA_DIR}/bin/conda config --set channel_priority strict \
    && ${CONDA_DIR}/bin/conda clean -afy

# -----------------------------------------------------------------------------
# CASCADE env: Protenix 1.0.4 + pipeline requirements
# -----------------------------------------------------------------------------
COPY requirements.txt /tmp/requirements.txt
RUN ${CONDA_DIR}/bin/conda create -n cascade python=3.11 pip -y -q \
    && /opt/conda/envs/cascade/bin/pip install --upgrade pip \
    && /opt/conda/envs/cascade/bin/pip install "protenix>=1.0.4,<2.0" \
    && /opt/conda/envs/cascade/bin/pip install -r /tmp/requirements.txt \
    && ${CONDA_DIR}/bin/conda clean -afy

# -----------------------------------------------------------------------------
# PXDesign env: clone upstream and run its install.sh, which creates the
# 'pxdesign' conda env with Protenix 0.5.0+pxd. Pinned CUDA 12.1 for wheel
# compatibility. This layer is the slow/large one; keep it separate from the
# code copy so application-code edits don't invalidate it.
# -----------------------------------------------------------------------------
ARG PXDESIGN_REF=main
RUN git clone --depth 1 --branch ${PXDESIGN_REF} \
        https://github.com/bytedance/PXDesign.git /opt/PXDesign \
    && cd /opt/PXDesign \
    && sed -i 's/python=3\.11/python=3.11 pip/g' install.sh \
    && bash -x install.sh --env pxdesign --pkg_manager conda --cuda-version 12.1 \
    && ${CONDA_DIR}/bin/conda clean -afy \
    && rm -rf /root/.cache

# -----------------------------------------------------------------------------
# Application code (last, so code edits don't bust the conda layers)
# -----------------------------------------------------------------------------
WORKDIR /workspace/CASCADE
COPY . /workspace/CASCADE

# Ensure entrypoint is executable and run-scoped dirs exist so the orchestrator
# never has to create its parent chain on first write.
RUN chmod +x /workspace/CASCADE/orchestrator/entrypoint.sh \
    && mkdir -p /workspace/CASCADE/outputs/runs \
                /workspace/CASCADE/logs

# Basic smoke test: the orchestrator module imports cleanly in the cascade env.
# Runs with cwd=scripts/ so the `from utils.protenix_eval import ...` resolves.
RUN cd /workspace/CASCADE/scripts && /opt/conda/envs/cascade/bin/python -c "\
import sys, importlib.util; \
sys.path.insert(0, '.'); \
spec=importlib.util.spec_from_file_location('eo','evolution_orchestrator.py'); \
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); \
print('orchestrator import OK, RUN_ID=', m.RUN_ID)"

ENTRYPOINT ["/usr/bin/tini", "--", "/workspace/CASCADE/orchestrator/entrypoint.sh"]
CMD []
