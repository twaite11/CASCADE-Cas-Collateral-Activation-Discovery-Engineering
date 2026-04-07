"""
Unit tests for 01_parse_and_annotate: DB init, FASTA/CSV load, HEPN detection, JSON gen.
No GPU required. Uses temp dirs for isolation.
"""
import json
import sqlite3
import pytest
from pathlib import Path

import importlib.util
import sys

# Load parse module
parse_module_path = Path(__file__).resolve().parent.parent / "scripts" / "01_parse_and_annotate.py"
spec = importlib.util.spec_from_file_location("parse_and_annotate", parse_module_path)
parse_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parse_mod)


def _patch_config(monkeypatch, data_dir, db_file, json_dir, meta_file):
    """Monkeypatch module constants for tests."""
    monkeypatch.setattr(parse_mod, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(parse_mod, "DB_FILE", str(db_file))
    monkeypatch.setattr(parse_mod, "JSON_OUT_DIR", str(json_dir))
    monkeypatch.setattr(parse_mod, "METADATA_OUT_FILE", str(meta_file))


class TestInitDb:
    """Test init_db."""

    def test_creates_schema(self, tmpdir, monkeypatch):
        db_file = tmpdir / "metadata" / "test.db"
        db_file.parent.mkdir(parents=True)
        _patch_config(monkeypatch, tmpdir, db_file, tmpdir / "jsons", tmpdir / "meta.json")
        conn = parse_mod.init_db()
        cur = conn.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variants'")
        assert cur.fetchone() is not None
        conn.close()


class TestLoadFilesToDb:
    """Test load_files_to_db."""

    def test_loads_fasta_and_csv(self, tmpdir, sample_fasta, sample_csv, monkeypatch):
        data_dir = Path(sample_fasta).parent
        db_file = tmpdir / "meta" / "test.db"
        db_file.parent.mkdir(parents=True)
        _patch_config(monkeypatch, data_dir, db_file, tmpdir / "jsons", tmpdir / "metadata.json")
        conn = parse_mod.init_db()
        count = parse_mod.load_files_to_db(conn)
        assert count >= 1
        cur = conn.cursor()
        cur.execute("SELECT sequence_id, crrna_repeat FROM variants WHERE sequence_id='test_cas13'")
        row = cur.fetchone()
        assert row is not None
        conn.close()

    def test_empty_dir_returns_zero(self, tmpdir, monkeypatch):
        empty = tmpdir / "empty"
        empty.mkdir()
        db_file = tmpdir / "test.db"
        _patch_config(monkeypatch, empty, db_file, tmpdir / "jsons", tmpdir / "meta.json")
        conn = parse_mod.init_db()
        count = parse_mod.load_files_to_db(conn)
        assert count == 0
        conn.close()


class TestIdentifyHepnDomains:
    """Test identify_hepn_domains."""

    def test_marks_hepn_success(self, tmpdir, sample_fasta, sample_csv, monkeypatch):
        data_dir = Path(sample_fasta).parent
        db_file = tmpdir / "meta" / "test.db"
        db_file.parent.mkdir(parents=True)
        _patch_config(monkeypatch, data_dir, db_file, tmpdir / "jsons", tmpdir / "meta.json")
        conn = parse_mod.init_db()
        parse_mod.load_files_to_db(conn)
        parse_mod.identify_hepn_domains(conn)
        cur = conn.cursor()
        cur.execute("SELECT status, hepn1_start, hepn2_start FROM variants WHERE sequence_id='test_cas13'")
        row = cur.fetchone()
        assert row[0] == "success"
        assert row[1] is not None and row[2] is not None
        conn.close()

    def test_marks_failed_when_few_motifs(self, tmpdir, monkeypatch):
        """Sequence with 0 or 1 HEPN motif should be marked failed."""
        data_dir = tmpdir / "data"
        data_dir.mkdir()
        (data_dir / "bad.fasta").write_text(">bad\nMGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")
        (data_dir / "bad_metadata.csv").write_text("sequence_id,repeat_domains,sra_accession,score\nbad,AAA,SRR1,0.5\n")
        db_file = tmpdir / "test.db"
        _patch_config(monkeypatch, data_dir, db_file, tmpdir / "jsons", tmpdir / "meta.json")
        conn = parse_mod.init_db()
        parse_mod.load_files_to_db(conn)
        parse_mod.identify_hepn_domains(conn)
        cur = conn.cursor()
        cur.execute("SELECT status FROM variants WHERE sequence_id='bad'")
        row = cur.fetchone()
        assert row[0] == "failed"
        conn.close()


class TestSelectHepnPair:
    """Test _select_hepn_pair with the improved separation/positional filtering."""

    def test_picks_correct_pair_over_first_last(self, tmpdir, monkeypatch):
        """Sequence with spurious R...H in NTD and C-term; real HEPN pair in the middle."""
        import re
        # Build a ~900 residue sequence with R...H hits at multiple positions.
        # Spurious hit at pos ~5 (NTD), real HEPN1 at ~200, real HEPN2 at ~500,
        # spurious hit at ~890 (C-term). First/last would pick 5 and 890.
        seq = "G" * 5 + "RAAAAH" + "G" * 189 + "RBBBBH" + "G" * 294 + "RCCCCCH" + "G" * 389 + "RDDDDDH" + "G" * 10
        motif = re.compile(r'R.{3,6}H')
        matches = list(motif.finditer(seq))
        assert len(matches) == 4, f"Expected 4 R...H hits, got {len(matches)}"

        pair = parse_mod._select_hepn_pair(seq, motif)
        assert pair is not None
        h1, h2 = pair
        sep = h2 - h1
        assert 150 <= sep <= 600, f"Separation {sep} outside expected range"
        assert h1 != matches[0].start(), "Should NOT pick the first (NTD) hit"
        assert h2 != matches[-1].start(), "Should NOT pick the last (C-term) hit"

    def test_returns_none_for_single_motif(self):
        import re
        seq = "G" * 100 + "RAAAAH" + "G" * 100
        motif = re.compile(r'R.{3,6}H')
        assert parse_mod._select_hepn_pair(seq, motif) is None

    def test_returns_none_for_no_motifs(self):
        import re
        seq = "G" * 500
        motif = re.compile(r'R.{3,6}H')
        assert parse_mod._select_hepn_pair(seq, motif) is None

    def test_fallback_when_no_positional_match(self):
        """Two motifs both in the first half, but with valid separation."""
        import re
        seq = "G" * 10 + "RAAAAH" + "G" * 200 + "RBBBBH" + "G" * 800
        motif = re.compile(r'R.{3,6}H')
        pair = parse_mod._select_hepn_pair(seq, motif)
        assert pair is not None, "Fallback path should still find a valid pair"

    def test_rejects_close_pairs(self):
        """Two motifs only 50 residues apart — should be rejected."""
        import re
        seq = "G" * 200 + "RAAAAH" + "G" * 44 + "RBBBBH" + "G" * 200
        motif = re.compile(r'R.{3,6}H')
        pair = parse_mod._select_hepn_pair(seq, motif)
        assert pair is None


class TestGenerateProtenixJsons:
    """Test generate_protenix_jsons."""

    def test_generates_json_and_metadata(self, tmpdir, sample_fasta, sample_csv, monkeypatch):
        data_dir = Path(sample_fasta).parent
        json_dir = tmpdir / "jsons"
        meta_file = tmpdir / "variant_domain_metadata.json"
        db_file = tmpdir / "test.db"
        _patch_config(monkeypatch, data_dir, db_file, json_dir, meta_file)
        conn = parse_mod.init_db()
        parse_mod.load_files_to_db(conn)
        parse_mod.identify_hepn_domains(conn)
        count = parse_mod.generate_protenix_jsons(conn)
        assert count >= 1
        assert meta_file.exists()
        meta = json.loads(meta_file.read_text())
        assert "test_cas13" in meta
        assert "domains" in meta["test_cas13"]
        assert "HEPN1" in meta["test_cas13"]["domains"]
        jsons = list(json_dir.glob("*.json"))
        assert len(jsons) >= 1
        conn.close()
