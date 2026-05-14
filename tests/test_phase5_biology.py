"""Phase 5 biology / pipeline correctness regression tests.

Covers (the actionable subset of) the original audit list:
  * B-13 SQLite-driven CANDIDATES (with legacy-dict fallback)
  * B-15 subtype-restricted DR library in Tier 3
"""
from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


# ===========================================================================
# B-15 subtype-restricted DR library
# ===========================================================================
class TestB15_SubtypeRestrictedDrs:
    def test_cas13a_only_returns_a_subfamily(self):
        import discover_crrna as dc
        drs = dc._drs_for_subtype("cas13a")
        assert drs, "cas13a should have known DRs"
        for name in drs:
            assert "Cas13a" in name, f"Tier 3 leaked non-a DR for cas13a: {name}"

    def test_cas13b_includes_bt_near_neighbour(self):
        import discover_crrna as dc
        drs = dc._drs_for_subtype("cas13b")
        # Cas13bt is a divergent VI-B and is treated as a near neighbour.
        assert any("Cas13b" in n and "bt" not in n.lower() for n in drs), \
            "cas13b screen should include native Cas13b DRs"
        assert any("Cas13bt" in n for n in drs), \
            "cas13b screen should also include Cas13bt (near-neighbour)"
        # ...but no Cas13a / Cas13d should sneak in.
        for n in drs:
            assert "Cas13a" not in n and "Cas13d" not in n

    def test_cas13d_isolated(self):
        import discover_crrna as dc
        drs = dc._drs_for_subtype("cas13d")
        for name in drs:
            assert "Cas13d" in name, f"Cas13d screen leaked: {name}"

    def test_unknown_subtype_returns_nothing_dangerous(self):
        import discover_crrna as dc
        drs = dc._drs_for_subtype("not-a-cas13")
        # Should return empty (or at least nothing from a different subtype):
        # the function defaults to its own name, which won't match any DR.
        assert all(n.startswith("not-a-cas13") for n in drs) or drs == {}

    def test_tier3_no_longer_cartesian(self):
        """The Tier 3 docstring is the source of truth for the new behaviour."""
        src = (SCRIPTS / "discover_crrna.py").read_text(encoding="utf-8")
        assert "subtype-restricted" in src
        # The old "every DR vs every candidate" log message must be gone.
        assert "Library size: {len(KNOWN_DR_LIBRARY)} DRs x {len(CANDIDATES)}" not in src


# ===========================================================================
# B-13 SQLite-driven CANDIDATES
# ===========================================================================
class TestB13_SqliteCandidates:
    def test_legacy_fallback_when_db_absent(self, monkeypatch, tmp_path):
        # Point DB_FILE at a non-existent location and re-evaluate the
        # loader -- the module-level CANDIDATES has already been built, so
        # we test the loader function directly.
        import discover_crrna as dc
        monkeypatch.setattr(dc, "DB_FILE", tmp_path / "no-such.db")
        assert dc._load_candidates_from_db() is None
        # The module-level CANDIDATES dict must still be non-empty thanks
        # to the legacy fallback.
        assert dc.CANDIDATES, "CANDIDATES should fall back to the legacy dict"
        assert all("subtype" in v for v in dc.CANDIDATES.values())

    def test_loader_reads_minimal_schema(self, monkeypatch, tmp_path):
        import discover_crrna as dc
        db = tmp_path / "cas13_variants.db"
        conn = sqlite3.connect(str(db))
        conn.execute(
            "CREATE TABLE variants (sequence_id TEXT PRIMARY KEY, "
            "subtype TEXT, elite INTEGER, contig_id TEXT, "
            "blastx_qstart INTEGER, blastx_qend INTEGER)"
        )
        conn.execute(
            "INSERT INTO variants VALUES "
            "('seq_a_elite', 'cas13a', 1, 'contig_1', 100, 200), "
            "('seq_b_not_elite', 'cas13b', 0, 'contig_2', 300, 400)"
        )
        conn.commit()
        conn.close()
        monkeypatch.setattr(dc, "DB_FILE", db)
        loaded = dc._load_candidates_from_db()
        assert loaded is not None
        # Only the elite row should appear.
        assert "seq_a_elite" in loaded
        assert "seq_b_not_elite" not in loaded
        assert loaded["seq_a_elite"]["subtype"] == "cas13a"
        assert loaded["seq_a_elite"]["contig"] == "contig_1"

    def test_loader_handles_missing_optional_columns(self, monkeypatch, tmp_path):
        """The schema may lack blastx_* columns on older databases."""
        import discover_crrna as dc
        db = tmp_path / "cas13_variants.db"
        conn = sqlite3.connect(str(db))
        conn.execute(
            "CREATE TABLE variants (sequence_id TEXT PRIMARY KEY, "
            "subtype TEXT)"
        )
        conn.execute("INSERT INTO variants VALUES ('seq_x', 'cas13d')")
        conn.commit()
        conn.close()
        monkeypatch.setattr(dc, "DB_FILE", db)
        loaded = dc._load_candidates_from_db()
        assert loaded is not None
        assert loaded["seq_x"]["subtype"] == "cas13d"
        assert loaded["seq_x"]["blastx_qstart"] == 0  # default

    def test_loader_returns_none_on_unusable_schema(self, monkeypatch, tmp_path):
        import discover_crrna as dc
        db = tmp_path / "cas13_variants.db"
        conn = sqlite3.connect(str(db))
        # Missing the required `subtype` column.
        conn.execute("CREATE TABLE variants (sequence_id TEXT PRIMARY KEY)")
        conn.execute("INSERT INTO variants VALUES ('seq_a')")
        conn.commit()
        conn.close()
        monkeypatch.setattr(dc, "DB_FILE", db)
        assert dc._load_candidates_from_db() is None
