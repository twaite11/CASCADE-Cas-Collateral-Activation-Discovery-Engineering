from __future__ import annotations

import logging
import os
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from .api import baselines as baselines_api
from .api import runs as runs_api
from .config import load_config
from .service import DashboardService
from .storage import SqliteVariantCatalogStore

log = logging.getLogger(__name__)

config = load_config()
catalog = SqliteVariantCatalogStore(config.sqlite_db_path)
service = DashboardService(config=config, catalog_store=catalog)

app = FastAPI(title="CASCADE Dashboard API", version="0.1.0")

# C-1 fix: never use the (allow_origins=*, allow_credentials=True) combo --
# modern browsers reject it outright AND any reverse proxy / edge that strips
# the wildcard turns this into a drive-by + CSRF surface.  Origins are read
# from CASCADE_ALLOWED_ORIGINS (comma-separated), defaulting to the two URLs
# the Vite dev server uses.  Setting CASCADE_ALLOWED_ORIGINS=* explicitly
# disables credentialed CORS so the upgrade path is clean.
_origins_env = os.environ.get(
    "CASCADE_ALLOWED_ORIGINS",
    "http://localhost:5173,http://127.0.0.1:5173",
).strip()
if _origins_env == "*":
    log.warning(
        "CASCADE_ALLOWED_ORIGINS=* — credentialed CORS is disabled; only "
        "anonymous read-only fetches from any origin will succeed."
    )
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=False,
        allow_methods=["GET", "HEAD", "OPTIONS"],
        allow_headers=["*"],
    )
else:
    _allowed_origins = [
        o.strip() for o in _origins_env.split(",") if o.strip()
    ]
    app.add_middleware(
        CORSMiddleware,
        allow_origins=_allowed_origins,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

app.include_router(baselines_api.router)
app.include_router(runs_api.router)
runs_api.register_dashboard_service(service)

_ui_dist = config.cascade_root / "dashboard_ui" / "dist"
_legacy_frontend = config.cascade_root / "dashboard_frontend"
frontend_dir = _ui_dist if _ui_dist.exists() else _legacy_frontend
if frontend_dir.exists():
    app.mount("/dashboard", StaticFiles(directory=frontend_dir, html=True), name="dashboard")


@app.get("/health")
def health() -> dict:
    return service.ping()


@app.get("/api/overview")
def overview() -> dict:
    return service.get_overview()


@app.get("/api/pipeline-health")
def pipeline_health() -> dict:
    return service.get_pipeline_health()


@app.get("/api/production")
def production() -> dict:
    return service.get_production_summary()


@app.get("/api/variants")
def variants(
    generation: int | None = Query(None),
    lineage: str | None = Query(None),
    optimized_only: bool = Query(False),
    elite_only: bool = Query(False),
    search: str | None = Query(None),
    min_fitness: float | None = Query(None),
    min_iptm: float | None = Query(None),
    min_af2_ig: float | None = Query(None),
    limit: int = Query(200, ge=1, le=2000),
    offset: int = Query(0, ge=0),
) -> dict:
    rows = service.load_variants()

    if generation is not None:
        rows = [r for r in rows if int(r.get("generation", -1) or -1) == generation]
    if optimized_only:
        rows = [r for r in rows if bool(r.get("optimized_switch", False))]
    if elite_only:
        rows = [r for r in rows if bool(r.get("is_elite", False))]
    if lineage:
        rows = [r for r in rows if str(r.get("baseline_id", "")) == lineage]
    if search:
        needle = search.lower()
        rows = [r for r in rows if needle in str(r.get("variant_id", "")).lower()]
    if min_fitness is not None:
        rows = [r for r in rows if float(r.get("fitness", -1e9) or -1e9) >= min_fitness]
    if min_iptm is not None:
        rows = [r for r in rows if float(r.get("iptm", 0.0) or 0.0) >= min_iptm]
    if min_af2_ig is not None:
        rows = [r for r in rows if float(r.get("af2_ig", 0.0) or 0.0) >= min_af2_ig]

    total = len(rows)
    page = rows[offset : offset + limit]
    return {"total": total, "offset": offset, "limit": limit, "rows": page}


@app.get("/api/optimized-switches")
def optimized_switches(limit: int = Query(500, ge=1, le=5000)) -> dict:
    rows = service.get_optimized_summary(limit=limit)
    return {"total": len(rows), "rows": rows}


@app.get("/api/variant/{variant_id}")
def variant_detail(variant_id: str) -> dict:
    row = service.get_variant(variant_id)
    if row is None:
        raise HTTPException(status_code=404, detail="Variant not found")
    return row


@app.get("/api/compare")
def compare_variants(left: str, right: str) -> dict:
    left_row = service.get_variant(left)
    right_row = service.get_variant(right)
    if left_row is None or right_row is None:
        raise HTTPException(status_code=404, detail="One or both variants not found")
    return {"left": left_row, "right": right_row}


@app.get("/api/structure-file")
def structure_file(path: str) -> FileResponse:
    # C-18 fix: previously used `Path.resolve()` which follows symlinks, then
    # checked containment by `root not in full_path.parents`.  A symlink
    # planted under cascade_root could point at `/etc/passwd` and pass the
    # check because the parents-of-the-target weren't compared.  Use
    # `os.path.commonpath` after `realpath` for true containment.
    try:
        candidate = os.path.realpath(str(config.cascade_root / path))
    except (OSError, ValueError):
        raise HTTPException(status_code=400, detail="Invalid path")
    root_real = os.path.realpath(str(config.cascade_root))
    try:
        common = os.path.commonpath([candidate, root_real])
    except ValueError:
        # Different drives on Windows -> definitely not contained.
        raise HTTPException(status_code=400, detail="Invalid path")
    if common != root_real:
        raise HTTPException(status_code=400, detail="Invalid path")
    full_path = Path(candidate)
    if not full_path.exists():
        raise HTTPException(status_code=404, detail="File not found")
    if full_path.suffix.lower() not in {".cif", ".pdb"}:
        raise HTTPException(status_code=400, detail="Unsupported structure type")
    return FileResponse(full_path)


@app.get("/", response_model=None)
def root_index() -> FileResponse | dict:
    index = frontend_dir / "index.html"
    if index.exists():
        return FileResponse(index)
    return {"message": "CASCADE dashboard API is running. Open /dashboard if frontend is present."}
