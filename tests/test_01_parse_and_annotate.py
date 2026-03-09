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
    """Test identify_hepn_domains with separation and positional filtering."""

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
        # HEPN2 should be well separated from HEPN1
        assert row[2] > row[1] + 100
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

    def test_marks_failed_when_motifs_too_close(self, tmpdir, monkeypatch):
        """Two HEPN motifs < 150 residues apart should be marked failed."""
        data_dir = tmpdir / "data"
        data_dir.mkdir()
        # Two R...H motifs only 100 residues apart
        close_seq = "MAAAAARAILXH" + "G" * 100 + "RVVVXH" + "G" * 40
        (data_dir / "close.fasta").write_text(f">close\n{close_seq}\n")
        (data_dir / "close_metadata.csv").write_text("sequence_id,repeat_domains,sra_accession,score\nclose,AAA,SRR1,0.5\n")
        db_file = tmpdir / "test.db"
        _patch_config(monkeypatch, data_dir, db_file, tmpdir / "jsons", tmpdir / "meta.json")
        conn = parse_mod.init_db()
        parse_mod.load_files_to_db(conn)
        parse_mod.identify_hepn_domains(conn)
        cur = conn.cursor()
        cur.execute("SELECT status FROM variants WHERE sequence_id='close'")
        row = cur.fetchone()
        assert row[0] == "failed"
        conn.close()

    def test_expanded_hepn_boundaries(self, tmpdir, sample_fasta, sample_csv, monkeypatch):
        """HEPN domains should use -60/+100 boundaries around the catalytic R."""
        data_dir = Path(sample_fasta).parent
        db_file = tmpdir / "meta" / "test.db"
        db_file.parent.mkdir(parents=True)
        _patch_config(monkeypatch, data_dir, db_file, tmpdir / "jsons", tmpdir / "meta.json")
        conn = parse_mod.init_db()
        parse_mod.load_files_to_db(conn)
        parse_mod.identify_hepn_domains(conn)
        cur = conn.cursor()
        cur.execute("SELECT hepn1_start, hepn1_end, hepn2_start, hepn2_end FROM variants WHERE sequence_id='test_cas13'")
        row = cur.fetchone()
        h1_start, h1_end, h2_start, h2_end = row
        # HEPN1 R is at pos 6: start = max(0, 6-60) = 0, end = 6+100 = 106
        assert h1_start == 0
        assert h1_end == 106
        # HEPN2 R is at pos 200: start = max(106, 200-60) = 140, end = 200+100 = 300
        assert h2_start == 140
        assert h2_end == 300
        conn.close()


class TestGenerateProtenixJsons:
    """Test generate_protenix_jsons (now generates OFF + ON JSONs)."""

    def test_generates_off_on_jsons_and_metadata(self, tmpdir, sample_fasta, sample_csv, monkeypatch):
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
        # Should have both OFF and ON JSONs
        off_jsons = list(json_dir.glob("*_OFF.json"))
        on_jsons = list(json_dir.glob("*_ON.json"))
        assert len(off_jsons) >= 1, "No OFF-state JSONs generated"
        assert len(on_jsons) >= 1, "No ON-state JSONs generated"
        conn.close()

    def test_off_json_has_no_target_rna(self, tmpdir, sample_fasta, sample_csv, monkeypatch):
        """OFF-state JSON should have protein + crRNA only (2 sequence entries)."""
        data_dir = Path(sample_fasta).parent
        json_dir = tmpdir / "jsons"
        meta_file = tmpdir / "variant_domain_metadata.json"
        db_file = tmpdir / "test.db"
        _patch_config(monkeypatch, data_dir, db_file, json_dir, meta_file)
        conn = parse_mod.init_db()
        parse_mod.load_files_to_db(conn)
        parse_mod.identify_hepn_domains(conn)
        parse_mod.generate_protenix_jsons(conn)
        off_jsons = list(json_dir.glob("*_OFF.json"))
        assert len(off_jsons) >= 1
        with open(off_jsons[0]) as f:
            payload = json.load(f)
        sequences = payload[0]["sequences"]
        assert len(sequences) == 2  # protein + crRNA only (no target RNA)
        assert "proteinChain" in sequences[0]
        assert "rnaSequence" in sequences[1]
        conn.close()

    def test_on_json_has_target_rna(self, tmpdir, sample_fasta, sample_csv, monkeypatch):
        """ON-state JSON should have protein + crRNA + target RNA (3 sequence entries)."""
        data_dir = Path(sample_fasta).parent
        json_dir = tmpdir / "jsons"
        meta_file = tmpdir / "variant_domain_metadata.json"
        db_file = tmpdir / "test.db"
        _patch_config(monkeypatch, data_dir, db_file, json_dir, meta_file)
        conn = parse_mod.init_db()
        parse_mod.load_files_to_db(conn)
        parse_mod.identify_hepn_domains(conn)
        parse_mod.generate_protenix_jsons(conn)
        on_jsons = list(json_dir.glob("*_ON.json"))
        assert len(on_jsons) >= 1
        with open(on_jsons[0]) as f:
            payload = json.load(f)
        sequences = payload[0]["sequences"]
        assert len(sequences) == 3  # protein + crRNA + target RNA
        assert "proteinChain" in sequences[0]
        assert "rnaSequence" in sequences[1]
        assert "rnaSequence" in sequences[2]
        conn.close()
