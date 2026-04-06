# CASCADE Dashboard Runbook

## Purpose

Read-only dashboard API + UI for CASCADE running on a VPS.

- Reads pipeline outputs from SQLite/JSON/JSONL/filesystem
- Exposes aggregated API endpoints
- Serves frontend dashboard with auto-refresh
- Does **not** mutate pipeline artifacts

## Start

```bash
cd /workspace/CASCADE
pip install -r requirements.txt
uvicorn dashboard_backend.main:app --host 0.0.0.0 --port 8000
```

Or:

```bash
./scripts/run_dashboard.sh
```

## Endpoints

- `GET /health`
- `GET /api/overview`
- `GET /api/pipeline-health`
- `GET /api/production`
- `GET /api/variants`
- `GET /api/optimized-switches`
- `GET /api/variant/{variant_id}`
- `GET /api/compare?left=<id>&right=<id>`
- `GET /api/structure-file?path=<relative-path>`
- `GET /dashboard`

## Read-only guarantees

- API layer performs no writes to:
  - `outputs/`
  - `metadata/`
  - pipeline scripts/data files
- Only filesystem reads and in-memory aggregation/caching are used.

## Environment variables

- `CASCADE_ROOT` (optional): override repository root used for data paths
- `DASHBOARD_HOST` (optional, `run_dashboard.sh`): bind host
- `DASHBOARD_PORT` (optional, `run_dashboard.sh`): bind port
- `DASH_MIN_IPTM` / `DASH_MIN_AF2_IG` / `DASH_MAX_ON_DISTANCE`: hybrid optimized-switch thresholds

## Troubleshooting

- `rl_dataset_missing` warning:
  - Ensure `outputs/rl_gym_data/rl_training_dataset.jsonl` exists.
- `sqlite_catalog_missing` warning:
  - Ensure `metadata/cas13_variants.db` exists (produced by ingest step).
- `rl_dataset_stale_over_6h`:
  - Pipeline likely not producing new records.
- Structure not loading in 3D panel:
  - Verify structure path exists and is under CASCADE root.
  - Supported formats: `.cif`, `.pdb`.
- Empty tables:
  - Confirm the pipeline has produced records in `rl_training_dataset.jsonl`.
  - Check `/health` and `/api/overview`.
