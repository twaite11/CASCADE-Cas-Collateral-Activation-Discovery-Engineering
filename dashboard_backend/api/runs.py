"""Runs REST + WebSocket API.

Wires the dashboard to the Vast.ai controller:

  * `GET /api/vast/offers`       — proxy to `vastai search offers`
  * `POST /api/runs`             — provision + launch in background
  * `GET /api/runs`              — list all runs
  * `GET /api/runs/{id}`         — detail + last-N log lines
  * `DELETE /api/runs/{id}`      — cancel + destroy instance
  * `WS /api/runs/{id}/logs`     — live log stream (xterm.js-compatible)

All state lives in SQLite + an in-process log hub; no Redis/Celery needed.
"""
from __future__ import annotations

import asyncio
import logging
import threading
from typing import Any

from fastapi import (
    APIRouter,
    Depends,
    HTTPException,
    Query,
    WebSocket,
    WebSocketDisconnect,
)
from pydantic import BaseModel, Field

from ..auth import require_api_key, require_api_key_ws
from ..config import DashboardConfig, load_config
from ..vast.controller import VastController
from ..vast.log_hub import LogHub, hub as _default_hub
from ..vast.promoter import promote_run
from ..vast.provisioner import VastCliError, VastProvisioner
from ..vast.runs_store import Run, RunsStore, RunStatus

log = logging.getLogger(__name__)

router = APIRouter(prefix="/api", tags=["runs"])


# ---------------------------------------------------------------------------
# Dependency graph (singletons bound to process lifetime)
# ---------------------------------------------------------------------------
# Option B (per-baseline fan-out) means three POST /api/runs can race the
# first-time lazy init of these singletons.  A simple module-wide lock
# serialises the check-then-create double-bookkeeping; subsequent calls
# hit the fast path under the lock with a sub-microsecond cost.
_state: dict[str, Any] = {}
_state_lock = threading.Lock()


def _config() -> DashboardConfig:
    return load_config()


def get_store(config: DashboardConfig = Depends(_config)) -> RunsStore:
    with _state_lock:
        store = _state.get("store")
        if store is None:
            store = RunsStore(config.runs_db_path)
            _state["store"] = store
    return store


def get_provisioner() -> VastProvisioner:
    with _state_lock:
        prov = _state.get("provisioner")
        if prov is None:
            prov = VastProvisioner()
            _state["provisioner"] = prov
    return prov


def get_hub() -> LogHub:
    return _default_hub


def get_controller(
    config: DashboardConfig = Depends(_config),
    store: RunsStore = Depends(get_store),
    provisioner: VastProvisioner = Depends(get_provisioner),
    hub: LogHub = Depends(get_hub),
) -> VastController:
    with _state_lock:
        ctrl = _state.get("controller")
        if ctrl is not None:
            return ctrl
        inputs_to_push = [
            config.cascade_root / "metadata",
            config.cascade_root / "jsons",
            config.cascade_root / "scripts",
            config.cascade_root / "outputs" / "phase1_screening",
            config.cascade_root / "outputs" / "validated_baseline_ids.txt",
        ]

        def _promote(run: Run, artifacts_dir, cascade_root) -> None:
            # Use the artifacts_dir's parent because controller passes the
            # run-specific directory; promote_run handles both shapes.
            promote_run(run.id, artifacts_dir, cascade_root)
            svc = _state.get("dashboard_service")
            if svc is not None and hasattr(svc, "invalidate"):
                svc.invalidate()

        ctrl = VastController(
            store=store,
            provisioner=provisioner,
            hub=hub,
            orchestrator_image=config.orchestrator_image,
            logs_dir=config.runs_logs_dir,
            artifacts_dir=config.runs_artifacts_dir,
            cascade_root=config.cascade_root,
            inputs_to_push=inputs_to_push,
            ssh_key_path=config.vast_ssh_key_path or None,
            promote_fn=_promote,
        )
        _state["controller"] = ctrl
    return ctrl


def register_dashboard_service(service: Any) -> None:
    """Let main.py hand us the long-lived DashboardService so the post-run
    promoter can invalidate its cache."""
    _state["dashboard_service"] = service


# ---------------------------------------------------------------------------
# Schemas
# ---------------------------------------------------------------------------
class LaunchRunRequest(BaseModel):
    baseline_ids: list[str] = Field(min_length=1)
    crrna_lookup_ids: list[str] | None = None
    offer_id: int
    max_generations: int = Field(12, ge=1, le=200)
    variants_per_gen: int = Field(5, ge=1, le=100)
    workers: int = Field(3, ge=1, le=16)
    label: str | None = None


class RunSummary(BaseModel):
    id: str
    label: str
    status: str
    created_at: float
    started_at: float | None
    finished_at: float | None
    cost_usd: float | None
    dph_usd: float | None
    baseline_ids: list[str]
    max_generations: int
    variants_per_gen: int
    instance_id: int | None
    ssh_host: str | None
    ssh_port: int | None
    error: str | None

    @classmethod
    def from_run(cls, run: Run) -> "RunSummary":
        d = run.to_dict()
        return cls(
            id=d["id"],
            label=d["label"],
            status=d["status"],
            created_at=d["created_at"],
            started_at=d["started_at"],
            finished_at=d["finished_at"],
            cost_usd=d["cost_usd"],
            dph_usd=d["dph_usd"],
            baseline_ids=d["baseline_ids"],
            max_generations=d["max_generations"],
            variants_per_gen=d["variants_per_gen"],
            instance_id=d["instance_id"],
            ssh_host=d["ssh_host"],
            ssh_port=d["ssh_port"],
            error=d["error"],
        )


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------
@router.get("/vast/offers")
async def list_offers(
    gpu_name: str = Query("A100_SXM4"),
    num_gpus: int = Query(1, ge=1, le=8),
    min_vram_gb: int = Query(40, ge=1),
    max_dph_total: float | None = Query(None),
    limit: int = Query(10, ge=1, le=50),
    provisioner: VastProvisioner = Depends(get_provisioner),
    _key: str = Depends(require_api_key),
) -> dict[str, Any]:
    # Vast.ai's search filter expects ``gpu_ram`` in GB (verified empirically:
    # the CLI returns ``VRAM=81.9`` for 80GB A100s and ``VRAM=41.0`` for 40GB
    # variants, and ``gpu_ram>=38912`` matches zero offers while
    # ``gpu_ram>=38`` matches both).  Previously we multiplied by 1024
    # assuming MB, which silently nuked every search.
    parts = [
        f"gpu_name={gpu_name}",
        f"num_gpus={num_gpus}",
        f"gpu_ram>={min_vram_gb}",
        "rentable=true",
        "verified=true",
    ]
    if max_dph_total is not None:
        parts.append(f"dph_total<={max_dph_total}")
    query = " ".join(parts)
    try:
        offers = await provisioner.search_offers(query=query, limit=limit)
    except VastCliError as exc:
        raise HTTPException(status_code=502, detail=str(exc))
    return {"query": query, "offers": [o.to_dict() for o in offers]}


@router.post("/runs", status_code=202)
async def launch_run(
    body: LaunchRunRequest,
    controller: VastController = Depends(get_controller),
    _key: str = Depends(require_api_key),
) -> RunSummary:
    if body.crrna_lookup_ids and len(body.crrna_lookup_ids) != len(body.baseline_ids):
        raise HTTPException(
            status_code=400,
            detail="crrna_lookup_ids must match baseline_ids length when provided",
        )
    try:
        run = await controller.launch(
            offer_id=body.offer_id,
            baseline_ids=body.baseline_ids,
            crrna_lookup_ids=body.crrna_lookup_ids,
            max_generations=body.max_generations,
            variants_per_gen=body.variants_per_gen,
            workers=body.workers,
            label=body.label,
        )
    except VastCliError as exc:
        raise HTTPException(status_code=502, detail=str(exc))
    return RunSummary.from_run(run)


@router.get("/runs")
def list_runs(
    status: list[str] | None = Query(None),
    limit: int = Query(200, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    store: RunsStore = Depends(get_store),
) -> dict[str, Any]:
    statuses = None
    if status:
        try:
            statuses = [RunStatus(s) for s in status]
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=f"invalid status: {exc}")
    runs = store.list(statuses=statuses, limit=limit, offset=offset)
    # C-4 fix: total used to be len(runs) which is just the current page;
    # ask the store for the real total count under the same status filter.
    total = store.count(statuses=statuses)
    return {
        "total": total,
        "offset": offset,
        "limit": limit,
        "rows": [RunSummary.from_run(r).model_dump() for r in runs],
    }


@router.get("/runs/{run_id}")
def get_run(
    run_id: str,
    tail_lines: int = Query(500, ge=0, le=5000),
    store: RunsStore = Depends(get_store),
    hub: LogHub = Depends(get_hub),
) -> dict[str, Any]:
    run = store.get(run_id)
    if not run:
        raise HTTPException(status_code=404, detail="run not found")
    snapshot = hub.snapshot(run_id)
    tail = snapshot[-tail_lines:] if tail_lines else []
    return {
        "run": RunSummary.from_run(run).model_dump(),
        "log_tail": tail,
    }


@router.delete("/runs/{run_id}", status_code=202)
async def delete_run(
    run_id: str,
    controller: VastController = Depends(get_controller),
    store: RunsStore = Depends(get_store),
    _key: str = Depends(require_api_key),
) -> dict[str, Any]:
    run = store.get(run_id)
    if not run:
        raise HTTPException(status_code=404, detail="run not found")
    await controller.cancel(run_id)
    return {"ok": True, "run_id": run_id}


@router.websocket("/runs/{run_id}/logs")
async def runs_ws_logs(websocket: WebSocket, run_id: str) -> None:
    # C-2 / C-25 fix: require an API key BEFORE accept() so unauthenticated
    # clients never receive any backend bytes; and use the shared dependency
    # graph (`get_store` / `get_hub`) instead of rebuilding the store with a
    # fresh load_config() call -- that bypassed the singleton cache and could
    # diverge from the one used by the HTTP routes.
    if not await require_api_key_ws(websocket):
        return
    await websocket.accept()
    hub = get_hub()
    store = get_store(load_config())
    run = store.get(run_id)
    if not run:
        await websocket.send_text("[server] run not found")
        await websocket.close(code=4004)
        return

    queue = await hub.subscribe(run_id)
    hb: asyncio.Task | None = None
    try:
        async def heartbeat() -> None:
            while True:
                await asyncio.sleep(20)
                try:
                    await websocket.send_text("\u001b[2m[heartbeat]\u001b[0m")
                except Exception:  # noqa: BLE001
                    return

        hb = asyncio.create_task(heartbeat())
        try:
            while True:
                line = await queue.get()
                await websocket.send_text(line)
        finally:
            # C-26 fix: cancel + await the heartbeat so we don't leak
            # "Task was destroyed but it is pending!" warnings.
            if hb is not None:
                hb.cancel()
                try:
                    await hb
                except (asyncio.CancelledError, Exception):  # noqa: BLE001
                    pass
    except WebSocketDisconnect:
        pass
    except Exception:  # noqa: BLE001
        log.exception("WS log stream error for %s", run_id)
    finally:
        await hub.unsubscribe(run_id, queue)
        try:
            await websocket.close()
        except Exception:  # noqa: BLE001
            pass
