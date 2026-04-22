# Orchestrator container

This directory packages the evolution orchestrator (`scripts/evolution_orchestrator.py`) as a self-contained GPU Docker image for parallel Vast.ai runs.

## What's inside the image

| Path                          | Contents                                                   |
|-------------------------------|------------------------------------------------------------|
| `/opt/conda/envs/cascade`     | Python 3.11 + Protenix 1.0.4 + CASCADE pipeline deps       |
| `/opt/conda/envs/pxdesign`    | PXDesign + Protenix 0.5.0+pxd (variant generation)         |
| `/workspace/CASCADE`          | This repository (code + small reference data)              |

Large runtime inputs — `jsons/`, `outputs/phase1_screening/`, `data/mined_hits/` — are **not** baked into the image. The controller rsyncs them onto each freshly provisioned VPS before kicking off the run.

## Build locally

```bash
docker build -f orchestrator/Dockerfile -t cascade-orchestrator:dev .
```

The build runs `PXDesign/install.sh`, which takes 10-20 minutes and produces a ~20 GB image. Increase `--pull` cache hits with BuildKit if iterating.

## Run manually (for debugging)

```bash
docker run --rm --gpus all \
  -e CASCADE_RUN_ID=local-test \
  -e CASCADE_BASELINE_IDS=baseline_A,baseline_B \
  -e CASCADE_MAX_GENERATIONS=3 \
  -v "$(pwd)/metadata:/workspace/CASCADE/metadata:ro" \
  -v "$(pwd)/jsons:/workspace/CASCADE/jsons:ro" \
  -v "$(pwd)/outputs/phase1_screening:/workspace/CASCADE/outputs/phase1_screening:ro" \
  -v "$(pwd)/outputs/runs:/workspace/CASCADE/outputs/runs" \
  cascade-orchestrator:dev
```

The entrypoint forwards any trailing `docker run ... <args>` through to the Python orchestrator, so CLI flags also work:

```bash
docker run --rm --gpus all cascade-orchestrator:dev \
  --run-id local-test --baseline-ids A,B --max-generations 3
```

## Run on Vast.ai

The dashboard controller provisions instances via the `vastai` CLI with `--onstart-cmd` set to the entrypoint invocation plus the chosen baselines. See `CONTAINER_DEPLOY.md` for the full recipe.

## Publishing

Tagged as `ghcr.io/<owner>/cascade-orchestrator`. CI workflow at `.github/workflows/build-orchestrator.yml` builds on main and on any `orchestrator-v*` tag.
