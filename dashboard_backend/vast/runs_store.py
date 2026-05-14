"""SQLite-backed store for dashboard-launched runs.

One row per evolution run. Status transitions:

    QUEUED -> PROVISIONING -> STARTING -> RUNNING -> SYNCING -> COMPLETED
                                                                 or -> FAILED
    any non-terminal state -> CANCELLED (user destroy)

Writes are protected by ``PRAGMA journal_mode=WAL`` so the background
asyncio runner and the FastAPI request handlers can operate concurrently.
"""
from __future__ import annotations

import json
import sqlite3
import time
import uuid
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Iterable


class RunStatus(str, Enum):
    QUEUED = "queued"
    PROVISIONING = "provisioning"
    STARTING = "starting"
    RUNNING = "running"
    SYNCING = "syncing"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


_TERMINAL_STATUSES = {RunStatus.COMPLETED, RunStatus.FAILED, RunStatus.CANCELLED}


@dataclass
class Run:
    id: str
    label: str
    status: RunStatus
    baseline_ids: list[str]
    crrna_lookup_ids: list[str]
    max_generations: int
    variants_per_gen: int
    workers: int
    offer_id: int | None
    instance_id: int | None = None
    ssh_host: str | None = None
    ssh_port: int | None = None
    ssh_user: str = "root"
    docker_image: str | None = None
    dph_usd: float | None = None
    created_at: float = field(default_factory=time.time)
    started_at: float | None = None
    finished_at: float | None = None
    cost_usd: float | None = None
    error: str | None = None
    logs_path: str | None = None
    artifacts_path: str | None = None

    @property
    def is_terminal(self) -> bool:
        return self.status in _TERMINAL_STATUSES

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "label": self.label,
            "status": self.status.value,
            "baseline_ids": list(self.baseline_ids),
            "crrna_lookup_ids": list(self.crrna_lookup_ids),
            "max_generations": self.max_generations,
            "variants_per_gen": self.variants_per_gen,
            "workers": self.workers,
            "offer_id": self.offer_id,
            "instance_id": self.instance_id,
            "ssh_host": self.ssh_host,
            "ssh_port": self.ssh_port,
            "ssh_user": self.ssh_user,
            "docker_image": self.docker_image,
            "dph_usd": self.dph_usd,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "cost_usd": self.cost_usd,
            "error": self.error,
            "logs_path": self.logs_path,
            "artifacts_path": self.artifacts_path,
        }


_SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    id TEXT PRIMARY KEY,
    label TEXT NOT NULL,
    status TEXT NOT NULL,
    baseline_ids TEXT NOT NULL,
    crrna_lookup_ids TEXT NOT NULL,
    max_generations INTEGER NOT NULL,
    variants_per_gen INTEGER NOT NULL,
    workers INTEGER NOT NULL,
    offer_id INTEGER,
    instance_id INTEGER,
    ssh_host TEXT,
    ssh_port INTEGER,
    ssh_user TEXT NOT NULL DEFAULT 'root',
    docker_image TEXT,
    dph_usd REAL,
    created_at REAL NOT NULL,
    started_at REAL,
    finished_at REAL,
    cost_usd REAL,
    error TEXT,
    logs_path TEXT,
    artifacts_path TEXT
);
CREATE INDEX IF NOT EXISTS idx_runs_status ON runs(status);
CREATE INDEX IF NOT EXISTS idx_runs_created ON runs(created_at);
"""


class RunsStore:
    """Thin SQLite wrapper; call sites use blocking sqlite3 off the event loop
    since reads/writes are tiny and WAL makes them cheap. If you need strict
    async, wrap calls in :func:`asyncio.to_thread`."""

    def __init__(self, db_path: Path) -> None:
        self.db_path = db_path
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._bootstrap()

    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path, timeout=5)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA journal_mode=WAL;")
        conn.execute("PRAGMA synchronous=NORMAL;")
        conn.execute("PRAGMA busy_timeout=5000;")
        return conn

    def _bootstrap(self) -> None:
        conn = self._connect()
        try:
            conn.executescript(_SCHEMA)
            conn.commit()
        finally:
            conn.close()

    @staticmethod
    def new_id() -> str:
        return uuid.uuid4().hex[:12]

    def create(self, run: Run) -> Run:
        conn = self._connect()
        try:
            conn.execute(
                """
                INSERT INTO runs (
                    id, label, status,
                    baseline_ids, crrna_lookup_ids,
                    max_generations, variants_per_gen, workers,
                    offer_id, instance_id, ssh_host, ssh_port, ssh_user,
                    docker_image, dph_usd,
                    created_at, started_at, finished_at,
                    cost_usd, error, logs_path, artifacts_path
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    run.id,
                    run.label,
                    run.status.value,
                    json.dumps(run.baseline_ids),
                    json.dumps(run.crrna_lookup_ids),
                    run.max_generations,
                    run.variants_per_gen,
                    run.workers,
                    run.offer_id,
                    run.instance_id,
                    run.ssh_host,
                    run.ssh_port,
                    run.ssh_user,
                    run.docker_image,
                    run.dph_usd,
                    run.created_at,
                    run.started_at,
                    run.finished_at,
                    run.cost_usd,
                    run.error,
                    run.logs_path,
                    run.artifacts_path,
                ),
            )
            conn.commit()
        finally:
            conn.close()
        return run

    def update(self, run_id: str, **fields: Any) -> None:
        if not fields:
            return
        if "status" in fields and isinstance(fields["status"], RunStatus):
            fields["status"] = fields["status"].value
        if "baseline_ids" in fields and isinstance(fields["baseline_ids"], list):
            fields["baseline_ids"] = json.dumps(fields["baseline_ids"])
        if "crrna_lookup_ids" in fields and isinstance(fields["crrna_lookup_ids"], list):
            fields["crrna_lookup_ids"] = json.dumps(fields["crrna_lookup_ids"])

        cols = ", ".join(f"{k} = ?" for k in fields)
        params = list(fields.values()) + [run_id]
        conn = self._connect()
        try:
            conn.execute(f"UPDATE runs SET {cols} WHERE id = ?", params)
            conn.commit()
        finally:
            conn.close()

    def get(self, run_id: str) -> Run | None:
        conn = self._connect()
        try:
            row = conn.execute("SELECT * FROM runs WHERE id = ?", (run_id,)).fetchone()
        finally:
            conn.close()
        return _row_to_run(row) if row else None

    def list(
        self,
        *,
        statuses: Iterable[RunStatus] | None = None,
        limit: int = 200,
        offset: int = 0,
    ) -> list[Run]:
        where = ""
        params: list[Any] = []
        if statuses:
            placeholders = ",".join("?" for _ in statuses)
            where = f"WHERE status IN ({placeholders})"
            params.extend(s.value for s in statuses)
        conn = self._connect()
        try:
            rows = conn.execute(
                f"SELECT * FROM runs {where} ORDER BY created_at DESC LIMIT ? OFFSET ?",
                params + [limit, offset],
            ).fetchall()
        finally:
            conn.close()
        return [_row_to_run(r) for r in rows]

    def count(
        self,
        *,
        statuses: Iterable[RunStatus] | None = None,
    ) -> int:
        """Total number of runs matching ``statuses`` (or all when None).

        C-4 fix: the runs API used `total = len(rows)` which only described
        the current page.  Pagination consumers (RunsPage, OptimizedSidebar)
        undercounted.  This method returns the real total.
        """
        where = ""
        params: list[Any] = []
        if statuses:
            placeholders = ",".join("?" for _ in statuses)
            where = f"WHERE status IN ({placeholders})"
            params.extend(s.value for s in statuses)
        conn = self._connect()
        try:
            row = conn.execute(
                f"SELECT COUNT(*) AS n FROM runs {where}", params
            ).fetchone()
        finally:
            conn.close()
        return int(row["n"] if row else 0)

    def active_runs(self) -> list[Run]:
        return self.list(
            statuses=(
                RunStatus.QUEUED,
                RunStatus.PROVISIONING,
                RunStatus.STARTING,
                RunStatus.RUNNING,
                RunStatus.SYNCING,
            ),
            limit=500,
        )


def _row_to_run(row: sqlite3.Row) -> Run:
    return Run(
        id=row["id"],
        label=row["label"],
        status=RunStatus(row["status"]),
        baseline_ids=json.loads(row["baseline_ids"] or "[]"),
        crrna_lookup_ids=json.loads(row["crrna_lookup_ids"] or "[]"),
        max_generations=row["max_generations"],
        variants_per_gen=row["variants_per_gen"],
        workers=row["workers"],
        offer_id=row["offer_id"],
        instance_id=row["instance_id"],
        ssh_host=row["ssh_host"],
        ssh_port=row["ssh_port"],
        ssh_user=row["ssh_user"] or "root",
        docker_image=row["docker_image"],
        dph_usd=row["dph_usd"],
        created_at=row["created_at"],
        started_at=row["started_at"],
        finished_at=row["finished_at"],
        cost_usd=row["cost_usd"],
        error=row["error"],
        logs_path=row["logs_path"],
        artifacts_path=row["artifacts_path"],
    )
