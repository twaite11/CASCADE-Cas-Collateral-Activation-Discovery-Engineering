"""Baselines API — aggregates every known enzyme baseline with its crRNA,
HEPN coords, validated-flag, and whether a Phase 1 structure exists. This
feeds the frontend baseline picker that drives ``POST /api/runs``."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from fastapi import APIRouter, Depends, HTTPException, Query

from ..config import DashboardConfig, load_config
from ..storage import SqliteVariantCatalogStore, VariantCatalogStore


router = APIRouter(prefix="/api", tags=["baselines"])


def _config() -> DashboardConfig:
    return load_config()


def _catalog(config: DashboardConfig = Depends(_config)) -> VariantCatalogStore:
    return SqliteVariantCatalogStore(config.sqlite_db_path)


@dataclass
class BaselineRecord:
    baseline_id: str
    subtype: str | None
    sequence_length: int | None
    hepn1_start: int | None
    hepn1_end: int | None
    hepn2_start: int | None
    hepn2_end: int | None
    crrna_repeat: str | None
    crrna_spacer: str | None
    crrna_lookup_id: str
    has_phase1_structure: bool
    validated: bool
    sra_accession: str | None
    score: float | None
    status: str | None
    reason: str | None
    base_json: str | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "baseline_id": self.baseline_id,
            "subtype": self.subtype,
            "sequence_length": self.sequence_length,
            "hepn1_start": self.hepn1_start,
            "hepn1_end": self.hepn1_end,
            "hepn2_start": self.hepn2_start,
            "hepn2_end": self.hepn2_end,
            "crrna_repeat": self.crrna_repeat,
            "crrna_spacer": self.crrna_spacer,
            "crrna_lookup_id": self.crrna_lookup_id,
            "has_phase1_structure": self.has_phase1_structure,
            "validated": self.validated,
            "sra_accession": self.sra_accession,
            "score": self.score,
            "status": self.status,
            "reason": self.reason,
            "base_json": self.base_json,
        }


def _safe_read_json(path: Path) -> Any | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _load_validated_ids(config: DashboardConfig) -> set[str]:
    path = config.cascade_root / "outputs" / "validated_baseline_ids.txt"
    if not path.exists():
        return set()
    return {
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }


def _extract_crrna(base_json_path: Path) -> tuple[str | None, str | None]:
    """Return (repeat, spacer) from a base Protenix input JSON, if available.

    The Phase 1 JSONs look like::

        [{"name": "...", "sequences": [
            {"proteinChain": {...}},
            {"rnaSequence": {"sequence": "<repeat+spacer>"}},
            {"rnaSequence": {"sequence": "<target>"}}
        ]}]

    We pull the first rnaSequence as the crRNA and heuristically split it on
    the canonical repeat prefix if present in the metadata.
    """
    data = _safe_read_json(base_json_path)
    if not data or not isinstance(data, list) or not data:
        return None, None
    seqs = data[0].get("sequences", [])
    for ent in seqs:
        rna = ent.get("rnaSequence") or ent.get("rna")
        if rna and rna.get("sequence"):
            crrna = rna["sequence"]
            return crrna, None
    return None, None


def _scan_phase1_structures(config: DashboardConfig) -> set[str]:
    """Return baseline_ids that have at least one Phase-1 predicted structure."""
    root = config.cascade_root / "outputs" / "phase1_screening"
    if not root.exists():
        return set()
    ids: set[str] = set()
    for child in root.iterdir():
        if not child.is_dir():
            continue
        name = child.name
        if name.endswith("_pred"):
            ids.add(name[: -len("_pred")])
        else:
            ids.add(name)
    return ids


def build_baselines(
    config: DashboardConfig, catalog: VariantCatalogStore
) -> list[BaselineRecord]:
    meta = _safe_read_json(config.domain_metadata_path) or {}
    if not isinstance(meta, dict):
        meta = {}

    validated_ids = _load_validated_ids(config)
    phase1_ids = _scan_phase1_structures(config)
    catalog_rows = {r.sequence_id: r for r in catalog.fetch_all_variants()}

    jsons_dir = config.cascade_root / "jsons"
    out: list[BaselineRecord] = []

    for baseline_id, info in meta.items():
        info = info or {}
        domains = info.get("domains") or {}
        h1 = domains.get("HEPN1") or {}
        h2 = domains.get("HEPN2") or {}

        base_json = jsons_dir / f"{baseline_id}.json"
        repeat = info.get("crRNA_repeat_used")
        spacer: str | None = None
        if base_json.exists():
            _cr, _sp = _extract_crrna(base_json)
            if _cr and not repeat:
                repeat = _cr
            # Derive spacer as trailing portion not matching repeat prefix
            if _cr and repeat and _cr.startswith(repeat):
                spacer = _cr[len(repeat):] or None
            elif _cr and repeat is None:
                spacer = _cr

        cat_row = catalog_rows.get(baseline_id)
        out.append(
            BaselineRecord(
                baseline_id=baseline_id,
                subtype=info.get("subtype"),
                sequence_length=info.get("sequence_length"),
                hepn1_start=h1.get("start"),
                hepn1_end=h1.get("end"),
                hepn2_start=h2.get("start"),
                hepn2_end=h2.get("end"),
                crrna_repeat=repeat,
                crrna_spacer=spacer,
                crrna_lookup_id=baseline_id,
                has_phase1_structure=baseline_id in phase1_ids,
                validated=baseline_id in validated_ids,
                sra_accession=getattr(cat_row, "sra_accession", None),
                score=getattr(cat_row, "score", None),
                status=getattr(cat_row, "status", None),
                reason=getattr(cat_row, "reason", None),
                base_json=str(base_json.relative_to(config.cascade_root))
                if base_json.exists()
                else None,
            )
        )

    out.sort(
        key=lambda r: (
            not r.validated,
            not r.has_phase1_structure,
            r.baseline_id,
        )
    )
    return out


@router.get("/baselines")
def list_baselines(
    search: str | None = Query(None),
    validated_only: bool = Query(False),
    with_structure_only: bool = Query(False),
    subtype: str | None = Query(None),
    limit: int = Query(500, ge=1, le=5000),
    offset: int = Query(0, ge=0),
    config: DashboardConfig = Depends(_config),
    catalog: VariantCatalogStore = Depends(_catalog),
) -> dict[str, Any]:
    rows = build_baselines(config, catalog)

    if search:
        s = search.lower()
        rows = [r for r in rows if s in r.baseline_id.lower()]
    if validated_only:
        rows = [r for r in rows if r.validated]
    if with_structure_only:
        rows = [r for r in rows if r.has_phase1_structure]
    if subtype:
        rows = [r for r in rows if (r.subtype or "") == subtype]

    total = len(rows)
    page = rows[offset : offset + limit]
    return {
        "total": total,
        "offset": offset,
        "limit": limit,
        "rows": [r.to_dict() for r in page],
    }


@router.get("/baselines/{baseline_id}")
def baseline_detail(
    baseline_id: str,
    config: DashboardConfig = Depends(_config),
    catalog: VariantCatalogStore = Depends(_catalog),
) -> dict[str, Any]:
    rows = build_baselines(config, catalog)
    for r in rows:
        if r.baseline_id == baseline_id:
            return r.to_dict()
    raise HTTPException(status_code=404, detail="Baseline not found")
