from __future__ import annotations

from .base import VariantCatalogRecord, VariantCatalogStore


class PostgresVariantCatalogStore(VariantCatalogStore):
    """
    Postgres adapter placeholder.

    The API layer already depends on the abstract VariantCatalogStore, so
    migration later only requires implementing this class and wiring DSN env.
    """

    def __init__(self, dsn: str) -> None:
        self.dsn = dsn

    def fetch_all_variants(self) -> list[VariantCatalogRecord]:
        raise NotImplementedError("Postgres adapter is not implemented yet.")

    def ping(self) -> dict[str, object]:
        return {
            "ok": False,
            "backend": "postgres",
            "reason": "not_implemented",
            "dsn": self.dsn,
        }
