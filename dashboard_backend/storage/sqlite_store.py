from __future__ import annotations

import sqlite3
from pathlib import Path

from .base import VariantCatalogRecord, VariantCatalogStore


class SqliteVariantCatalogStore(VariantCatalogStore):
    def __init__(self, db_path: Path) -> None:
        self.db_path = db_path

    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path, timeout=5)
        conn.row_factory = sqlite3.Row
        # WAL significantly improves read/write coexistence on a VPS
        conn.execute("PRAGMA journal_mode=WAL;")
        conn.execute("PRAGMA synchronous=NORMAL;")
        conn.execute("PRAGMA busy_timeout=5000;")
        return conn

    def fetch_all_variants(self) -> list[VariantCatalogRecord]:
        if not self.db_path.exists():
            return []

        with self._connect() as conn:
            rows = conn.execute(
                """
                SELECT
                    sequence_id,
                    sra_accession,
                    score,
                    hepn1_start,
                    hepn1_end,
                    hepn2_start,
                    hepn2_end,
                    status,
                    reason
                FROM variants
                """
            ).fetchall()

        return [
            VariantCatalogRecord(
                sequence_id=row["sequence_id"],
                sra_accession=row["sra_accession"],
                score=row["score"],
                hepn1_start=row["hepn1_start"],
                hepn1_end=row["hepn1_end"],
                hepn2_start=row["hepn2_start"],
                hepn2_end=row["hepn2_end"],
                status=row["status"],
                reason=row["reason"],
            )
            for row in rows
        ]

    def ping(self) -> dict[str, object]:
        if not self.db_path.exists():
            return {"ok": False, "reason": "db_missing", "path": str(self.db_path)}

        try:
            with self._connect() as conn:
                conn.execute("SELECT 1").fetchone()
            return {"ok": True, "backend": "sqlite", "path": str(self.db_path)}
        except Exception as exc:  # noqa: BLE001
            return {"ok": False, "backend": "sqlite", "error": str(exc)}
