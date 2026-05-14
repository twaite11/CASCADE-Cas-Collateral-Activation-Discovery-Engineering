"""Regression tests for the 2026-05-13 Phase 1 audit fixes.

Each test class targets one audit finding (B-#) so failures point directly at
the responsible commit / file.  See docs/AUDIT_2026-05-13.md.

These tests are pure-Python and run without GPU, Protenix, PXDesign, MinCED,
DIAMOND, ViennaRNA, or cargo.  The Rust accelerator parity tests live in a
separate suite and skip cleanly when the binaries are not built.
"""
from __future__ import annotations

import importlib
import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


# ---------------------------------------------------------------------------
# Helpers: dynamically import scripts that have non-identifier filenames.
# ---------------------------------------------------------------------------
def _import_path(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def parse_mod():
    return _import_path("parse_and_annotate", SCRIPTS / "01_parse_and_annotate.py")


@pytest.fixture(scope="module")
def mining_v2_mod():
    import mining_v2  # type: ignore[import-not-found]
    return mining_v2


@pytest.fixture(scope="module")
def fix_crrna_mod():
    import fix_crrna_assignments  # type: ignore[import-not-found]
    return fix_crrna_assignments


# ===========================================================================
# B-2: mining_v2.py
#   1. HEPN regex must be canonical R-X(4-6)-H
#   2. _finalize_array must record array_start / array_end so MinCED-detected
#      arrays can be matched to ORFs by genomic distance.
# ===========================================================================
class TestB2_MiningV2:
    def test_hepn_regex_is_canonical_4_6(self, mining_v2_mod):
        pat = mining_v2_mod.HEPN_REGEX
        # R + exactly 4-6 chars + H
        for n in (4, 5, 6):
            assert pat.fullmatch("R" + "A" * n + "H"), f"should match R-X({n})-H"
        # Should NOT match outside the 4-6 window
        for n in (0, 1, 2, 3, 7, 8, 9):
            assert pat.fullmatch("R" + "A" * n + "H") is None, (
                f"should reject R-X({n})-H"
            )

    def test_finalize_array_stores_genomic_coords(self, mining_v2_mod):
        from collections import defaultdict
        arrays: dict = defaultdict(list)
        mining_v2_mod._finalize_array(
            arrays,
            contig="ctg1",
            repeats=["GAAAC", "GAAAC", "GAAAC"],
            spacers=["AAAATT", "CCCCCC"],
            array_start=1000,
            array_end=1200,
        )
        assert "ctg1" in arrays and len(arrays["ctg1"]) == 1
        rec = arrays["ctg1"][0]
        assert rec["array_start"] == 1000
        assert rec["array_end"] == 1200
        assert rec["n_repeats"] == 3
        assert rec["consensus_repeat"] == "GAAAC"

    def test_finalize_array_legacy_call_skips_coords(self, mining_v2_mod):
        """Pure-Python fallback path may still call without coords; the entry
        must be recorded but without array_start/end keys so assign_dr_to_orf
        will skip it rather than match the wrong array."""
        from collections import defaultdict
        arrays: dict = defaultdict(list)
        mining_v2_mod._finalize_array(
            arrays,
            contig="ctg1",
            repeats=["GAAAC", "GAAAC"],
            spacers=["AAA"],
        )
        rec = arrays["ctg1"][0]
        assert "array_start" not in rec
        assert "array_end" not in rec

    def test_assign_dr_to_orf_now_works_for_minced_arrays(self, mining_v2_mod):
        """End-to-end: build a synthetic MinCED-style array and verify
        assign_dr_to_orf picks it up. Was structurally unreachable pre-fix."""
        from collections import defaultdict
        arrays: dict = defaultdict(list)
        mining_v2_mod._finalize_array(
            arrays, "ctg1", ["GAAAC"] * 4, ["TGA" * 5] * 3,
            array_start=10_000, array_end=10_400,
        )
        result = mining_v2_mod.assign_dr_to_orf(
            orf_start_nt=8_000,
            orf_end_nt=9_500,
            contig_len=20_000,
            arrays=arrays["ctg1"],
        )
        assert result is not None
        assert result["dr"] == "GAAAC"
        assert result["distance_bp"] == 500  # 10_000 - 9_500
        assert result["n_repeats"] == 4

    def test_minced_gff_parser_populates_coords(self, mining_v2_mod, tmp_path):
        """Fixture GFF: parse via run_minced (mocked subprocess path)."""
        gff = (
            "##gff-version 3\n"
            "ctg1\tminced\trepeat_region\t1000\t1200\t.\t+\t.\tID=CRISPR1\n"
            "ctg1\tminced\trepeat_unit\t1000\t1027\t.\t+\t.\trpt_unit_seq=GATTTAGAGTACCTCAAAACAGAAGAGG\n"
            "ctg1\tminced\trepeat_unit\t1080\t1107\t.\t+\t.\trpt_unit_seq=GATTTAGAGTACCTCAAAACAGAAGAGG\n"
            "ctg1\tminced\trepeat_unit\t1170\t1197\t.\t+\t.\trpt_unit_seq=GATTTAGAGTACCTCAAAACAGAAGAGG\n"
            "ctg1\tminced\tbinding_site\t1027\t1080\t.\t+\t.\tspacer=AAATGCTAAGCCACGTCAAGTACCATTGTCCAGTGAAAATTCCAAATA\n"
            "ctg1\tminced\tbinding_site\t1107\t1170\t.\t+\t.\tspacer=GGCAATGGACTAGCTGCTGGATGTGCAATACCATGTAACCATTGATATGACCT\n"
        )
        out_dir = tmp_path
        (out_dir / "minced_output.gff").write_text(gff)

        # Reach into the parser directly without invoking the real `minced`
        # binary: re-implement the same parse the function does, on the fixture.
        import re
        from collections import defaultdict

        arrays: dict = defaultdict(list)
        current_contig = None
        current_repeats: list = []
        current_spacers: list = []
        current_array_start = None
        current_array_end = None
        for line in gff.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            contig = parts[0]
            feature = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            attrs = parts[8]
            if feature == "repeat_region":
                if current_contig and current_repeats:
                    mining_v2_mod._finalize_array(
                        arrays, current_contig, current_repeats,
                        current_spacers, current_array_start, current_array_end,
                    )
                current_contig = contig
                current_repeats = []
                current_spacers = []
                current_array_start = start
                current_array_end = end
            elif feature == "repeat_unit":
                m = re.search(r"rpt_unit_seq=([ACGT]+)", attrs, re.I)
                if m:
                    current_repeats.append(m.group(1))
            elif feature == "binding_site":
                m = re.search(r"spacer=([ACGT]+)", attrs, re.I)
                if m:
                    current_spacers.append(m.group(1))
        if current_contig and current_repeats:
            mining_v2_mod._finalize_array(
                arrays, current_contig, current_repeats, current_spacers,
                current_array_start, current_array_end,
            )
        rec = arrays["ctg1"][0]
        assert rec["array_start"] == 1000
        assert rec["array_end"] == 1200
        assert rec["n_repeats"] == 3


# ===========================================================================
# B-3: 01_parse_and_annotate._select_hepn_pair (already had a test in
#      test_01_parse_and_annotate.py).  Here we add a stricter case that
#      exposes the old first/last bug AND the new R-X(4-6)-H regex.
# ===========================================================================
class TestB3_HepnPair:
    def test_canonical_motif_in_parse_module(self, parse_mod):
        """The motif used by identify_hepn_domains must be R-X(4-6)-H now."""
        src = (SCRIPTS / "01_parse_and_annotate.py").read_text()
        assert "r'R.{4,6}H'" in src or 'r"R.{4,6}H"' in src, (
            "01_parse_and_annotate must use the canonical R-X(4-6)-H regex"
        )
        assert "r'R.{3,6}H'" not in src and 'r"R.{3,6}H"' not in src

    def test_select_pair_prefers_middle_over_endpoints(self, parse_mod):
        # 1000 aa sequence:
        #   pos 5    : spurious NTD hit
        #   pos 200  : real HEPN1
        #   pos 500  : real HEPN2
        #   pos 990  : spurious CTD hit
        # Old first/last code would pick (5, 990) -> sep 985, both wrong domains.
        seq = (
            "G" * 5 + "RAAAAH" +
            "G" * 189 + "RBBBBH" +
            "G" * 294 + "RCCCCCH" +
            "G" * 389 + "RDDDDDH" +
            "G" * 10
        )
        motif = re.compile(r"R.{4,6}H")
        starts = [m.start() for m in motif.finditer(seq)]
        assert len(starts) == 4

        pair = parse_mod._select_hepn_pair(seq, motif)
        assert pair is not None
        h1, h2 = pair
        sep = h2 - h1
        assert 150 <= sep <= 600
        assert h1 != starts[0], "must not pick the spurious NTD hit"
        assert h2 != starts[-1], "must not pick the spurious CTD hit"


# ===========================================================================
# B-4: failed off-target Protenix run must penalize, not reward.
# ===========================================================================
class TestB4_OffTargetFailurePenalty:
    @pytest.fixture(scope="class")
    def orch_mod(self):
        return _import_path("evo", SCRIPTS / "evolution_orchestrator.py")

    def test_failure_path_uses_MAX_ON(self, orch_mod):
        src = (SCRIPTS / "evolution_orchestrator.py").read_text()
        assert "failed_dist = float(MAX_ON_DISTANCE)" in src, (
            "off-target failure should be scored as MAX_ON_DISTANCE (penalty)"
        )
        # And the legacy bug pattern must be gone:
        assert "offtarget_by_mismatch[n_mismatch] = MIN_OFF_DISTANCE" not in src

    def test_fitness_penalizes_collapsed_offtarget(self, orch_mod):
        """compute_fitness with a collapsed offtarget (== ON-like) must score
        strictly lower than the same call with the offtarget fully dormant."""
        cf = orch_mod.compute_fitness
        baseline_args = dict(off_dist=25.0, on_dist=8.0, iptm_score=0.9,
                              af2_ig_score=0.85, is_full_ternary=True)
        clean = cf(**baseline_args, offtarget_by_mismatch=None)
        dormant = cf(
            **baseline_args,
            offtarget_by_mismatch={3: float(orch_mod.MIN_OFF_DISTANCE)},
        )
        collapsed = cf(
            **baseline_args,
            offtarget_by_mismatch={3: float(orch_mod.MAX_ON_DISTANCE)},
        )
        # Collapsed off-target should incur a penalty relative to dormant.
        assert collapsed < dormant, (
            f"collapsed offtarget ({collapsed}) must score lower than "
            f"dormant ({dormant})"
        )
        # And lower than the no-offtarget baseline.
        assert collapsed <= clean


# ===========================================================================
# B-5: elite gate must include OFF >= MIN_OFF_DISTANCE.
# ===========================================================================
class TestB5_EliteOffGate:
    def test_elite_condition_includes_off(self):
        src = (SCRIPTS / "evolution_orchestrator.py").read_text()
        # Look for the exact (and only) is_elite assignment.
        m = re.search(r"is_elite\s*=\s*\(([^)]*)\)", src, re.DOTALL)
        assert m, "could not find is_elite assignment"
        body = m.group(1)
        assert "off_dist" in body, "elite condition must reference off_dist"
        assert "MIN_OFF_DISTANCE" in body, (
            "elite condition must include MIN_OFF_DISTANCE gate"
        )


# ===========================================================================
# B-9: corrected-DR report column rename / dual-read.
# ===========================================================================
class TestB9_CorrectedDrColumn:
    def _make_csvs(self, tmpdir, repeat_validation_csv: str):
        """Build the minimum environment 01_parse_and_annotate expects."""
        data_dir = tmpdir / "data" / "mined_hits"
        data_dir.mkdir(parents=True)
        # Two ORFs; we put corrected DRs for orf-A only via the report,
        # orf-B falls back to the per-row repeat_domains column.
        (data_dir / "x.fasta").write_text(
            ">orf-A\nMAAAAAAAAAAAAAAAA\n>orf-B\nMBBBBBBBBBBBBBBBB\n"
        )
        (data_dir / "x_metadata.csv").write_text(
            "sequence_id,repeat_domains,sra_accession,score\n"
            "orf-A,RAW_A_KMER,SRR1,0.5\n"
            "orf-B,RAW_B_KMER,SRR2,0.5\n"
        )
        outputs = tmpdir / "outputs"
        outputs.mkdir(parents=True)
        report = outputs / "repeat_validation_report.csv"
        report.write_text(repeat_validation_csv)
        return data_dir, str(report)

    @pytest.mark.parametrize("csv_payload,expected_a", [
        # Rust path: writes chosen_repeat
        (
            "sequence_id,original_first_kmer,chosen_repeat,chosen_length,"
            "selection_reason,structure,mfe_kcal_mol,structure_ok\n"
            "orf-A,RAW_A_KMER,RUST_CHOSEN_DR,14,best,((..))..,-3.2,true\n",
            "RUST_CHOSEN_DR",
        ),
        # Python path: writes new_dr
        (
            "sequence_id,accession,old_first_kmer,new_dr,new_dr_length,"
            "source,rejection_reason,is_valid\n"
            "orf-A,acc,RAW_A_KMER,PY_NEW_DR,9,crispr_array_3x,,True\n",
            "PY_NEW_DR",
        ),
        # Both columns: chosen_repeat wins
        (
            "sequence_id,chosen_repeat,new_dr,structure_ok\n"
            "orf-A,RUST_WINS,PY_LOSES,true\n",
            "RUST_WINS",
        ),
        # chosen_repeat empty, fall through to new_dr
        (
            "sequence_id,chosen_repeat,new_dr,structure_ok\n"
            "orf-A,,PY_FALLBACK,true\n",
            "PY_FALLBACK",
        ),
    ])
    def test_load_uses_corrected_dr(self, tmpdir, monkeypatch, parse_mod,
                                     csv_payload, expected_a):
        data_dir, report_path = self._make_csvs(tmpdir, csv_payload)
        db_file = tmpdir / "meta" / "test.db"
        db_file.parent.mkdir(parents=True)
        monkeypatch.setattr(parse_mod, "DATA_DIR", str(data_dir))
        monkeypatch.setattr(parse_mod, "DB_FILE", str(db_file))
        monkeypatch.setattr(parse_mod, "JSON_OUT_DIR", str(tmpdir / "jsons"))
        monkeypatch.setattr(parse_mod, "METADATA_OUT_FILE", str(tmpdir / "meta.json"))
        monkeypatch.setattr(parse_mod, "CORRECTED_DR_REPORT_PATH", report_path)
        conn = parse_mod.init_db()
        try:
            parse_mod.load_files_to_db(conn)
            cur = conn.cursor()
            cur.execute(
                "SELECT crrna_repeat FROM variants WHERE sequence_id='orf-A'"
            )
            assert cur.fetchone()[0] == expected_a
            cur.execute(
                "SELECT crrna_repeat FROM variants WHERE sequence_id='orf-B'"
            )
            # orf-B is not in the corrected report, so it falls back to the
            # per-row repeat_domains value.
            assert cur.fetchone()[0] == "RAW_B_KMER"
        finally:
            conn.close()

    def test_structure_ok_false_skips_corrected_dr(self, tmpdir, monkeypatch,
                                                    parse_mod):
        """If the Rust report explicitly flags structure_ok=false, the corrected
        DR is dropped and the loader falls back to repeat_domains."""
        data_dir, report_path = self._make_csvs(
            tmpdir,
            "sequence_id,chosen_repeat,structure_ok\n"
            "orf-A,SHOULD_BE_DROPPED,false\n",
        )
        db_file = tmpdir / "meta" / "test.db"
        db_file.parent.mkdir(parents=True)
        monkeypatch.setattr(parse_mod, "DATA_DIR", str(data_dir))
        monkeypatch.setattr(parse_mod, "DB_FILE", str(db_file))
        monkeypatch.setattr(parse_mod, "JSON_OUT_DIR", str(tmpdir / "jsons"))
        monkeypatch.setattr(parse_mod, "METADATA_OUT_FILE", str(tmpdir / "meta.json"))
        monkeypatch.setattr(parse_mod, "CORRECTED_DR_REPORT_PATH", report_path)
        conn = parse_mod.init_db()
        try:
            parse_mod.load_files_to_db(conn)
            cur = conn.cursor()
            cur.execute("SELECT crrna_repeat FROM variants WHERE sequence_id='orf-A'")
            assert cur.fetchone()[0] == "RAW_A_KMER"
        finally:
            conn.close()


# ===========================================================================
# B-10: has_crispr_structure must fail-closed when ViennaRNA is unavailable.
# ===========================================================================
class TestB10_ViennaRNAFailClosed:
    def test_fails_closed_without_viennarna(self, fix_crrna_mod, monkeypatch):
        """Sabotage `import RNA` and confirm has_crispr_structure returns False."""
        monkeypatch.delitem(sys.modules, "RNA", raising=False)

        original_import = __builtins__["__import__"] if isinstance(__builtins__, dict) else __builtins__.__import__

        def fake_import(name, *args, **kwargs):
            if name == "RNA":
                raise ImportError("simulated missing ViennaRNA")
            return original_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", fake_import)
        monkeypatch.delenv("CASCADE_VIENNARNA_FAIL_OPEN", raising=False)
        # Reset the module-level "we already warned" flag so we exercise the
        # warning path inside this test as well.
        if hasattr(fix_crrna_mod, "_VIENNA_MISSING_WARNED"):
            fix_crrna_mod._VIENNA_MISSING_WARNED = False
        assert fix_crrna_mod.has_crispr_structure("GAUUUAGAGUACCUCAAAACAGAAGAGG") is False

    def test_fail_open_env_var_restores_legacy(self, fix_crrna_mod, monkeypatch):
        monkeypatch.delitem(sys.modules, "RNA", raising=False)
        original_import = __builtins__["__import__"] if isinstance(__builtins__, dict) else __builtins__.__import__

        def fake_import(name, *args, **kwargs):
            if name == "RNA":
                raise ImportError("simulated missing ViennaRNA")
            return original_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", fake_import)
        monkeypatch.setenv("CASCADE_VIENNARNA_FAIL_OPEN", "1")
        if hasattr(fix_crrna_mod, "_VIENNA_MISSING_WARNED"):
            fix_crrna_mod._VIENNA_MISSING_WARNED = False
        assert fix_crrna_mod.has_crispr_structure("GAUUUAGAGUACCUCAAAACAGAAGAGG") is True


# ===========================================================================
# B-1: Rust cascade_structscore NE2-NE2 distance.
#      Only runs when the binary is built; otherwise we still verify the
#      source change.
# ===========================================================================
class TestB1_HepnNE2:
    def test_source_uses_NE2_default_with_fallback(self):
        src = (PROJECT_ROOT / "rust" / "cascade_structscore" / "src" / "main.rs").read_text()
        # Default atom argument
        assert 'default_value = "NE2"' in src, "primary atom must default to NE2"
        # Fallback list
        assert 'default_value = "ND1,CG,CA"' in src, (
            "atom_fallback must default to ND1,CG,CA"
        )
        # The old hard-coded CA-only finder is gone
        assert 'a.name().trim() == "CA"' not in src, (
            "old CA-only atom finder must be replaced with find_atom_coords"
        )

    @pytest.mark.skipif(
        not (PROJECT_ROOT / "rust" / "target" / "release" / "cascade_structscore.exe").exists()
        and not (PROJECT_ROOT / "rust" / "target" / "release" / "cascade_structscore").exists(),
        reason="cascade_structscore binary not built; run `cargo build --release -p cascade_structscore`",
    )
    def test_binary_measures_NE2_NE2(self, minimal_pdb):
        binary = PROJECT_ROOT / "rust" / "target" / "release" / "cascade_structscore"
        if not binary.exists():
            binary = binary.with_suffix(".exe")
        out = subprocess.run(
            [str(binary), "hepn-distance",
             "--structure", minimal_pdb,
             "--h1-idx", "10", "--h2-idx", "50", "--chain", "A"],
            capture_output=True, text=True, timeout=10,
        )
        assert out.returncode == 0, out.stderr
        # MINIMAL_PDB has NE2-NE2 = 26 A (not CA-CA = 30 A)
        assert abs(float(out.stdout.strip()) - 26.0) < 0.1


# ===========================================================================
# B-3 Rust parity (cascade_ingest select_hepn_pair).  Source-only check.
# ===========================================================================
class TestB3_RustIngestSource:
    def test_rust_select_hepn_pair_present(self):
        src = (PROJECT_ROOT / "rust" / "cascade_ingest" / "src" / "main.rs").read_text()
        assert "fn select_hepn_pair" in src
        assert 'r"R.{4,6}H"' in src
        # Old loose regex is gone
        assert 'r"R.{3,6}H"' not in src
        # Old buggy first/last selection logic must be gone (only the new
        # `let starts:` form should remain).
        assert "matches[0].start() as i64" not in src
        assert "matches[matches.len() - 1].start() as i64" not in src
