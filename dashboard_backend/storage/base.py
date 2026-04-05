from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class VariantCatalogRecord:
    sequence_id: str
    sra_accession: str | None
    score: float | None
    hepn1_start: int | None
    hepn1_end: int | None
    hepn2_start: int | None
    hepn2_end: int | None
    status: str | None
    reason: str | None


class VariantCatalogStore:
    def fetch_all_variants(self) -> list[VariantCatalogRecord]:
        raise NotImplementedError

    def ping(self) -> dict[str, Any]:
        raise NotImplementedError
