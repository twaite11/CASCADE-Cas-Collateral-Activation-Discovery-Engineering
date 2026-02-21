"""
Unit tests for protenix_eval: JSON generation, mismatch sequences. No GPU/subprocess.
"""
import json
import pytest
from pathlib import Path

from utils.protenix_eval import (
    TARGET_REGION,
    DUMMY_SPACER_RNA,
    DUMMY_TARGET_RNA,
    generate_offtarget_sequences,
    generate_mismatch_sequences,
    generate_offtarget_json,
    generate_evaluation_jsons,
)


class TestGenerateMismatchSequences:
    """Test generate_mismatch_sequences."""

    def test_returns_correct_count(self):
        results = generate_mismatch_sequences(TARGET_REGION, mismatch_counts=(1, 2, 3), num_per_count=1, seed=42)
        assert len(results) == 3  # 1 + 2 + 3 mismatch
        assert all(isinstance(r, tuple) and len(r) == 2 for r in results)
        for rna, n in results:
            assert n in (1, 2, 3)
            assert rna.startswith("AAAAAA") and rna.endswith("AAAAAA")
            assert len(rna) == len(DUMMY_TARGET_RNA)

    def test_seed_reproducibility(self):
        a = generate_mismatch_sequences(TARGET_REGION, seed=123)
        b = generate_mismatch_sequences(TARGET_REGION, seed=123)
        assert a == b

    def test_exact_mismatch_count(self):
        """Each sequence should differ from target by exactly n positions."""
        target = "ACGUACGUACGUACGUCAGUCGAC"
        for rna, n_mm in generate_mismatch_sequences(target, (1, 2), num_per_count=1, seed=0):
            core = rna[6:-6]  # strip AAAAAA flank
            diffs = sum(1 for a, b in zip(core, target) if a != b)
            assert diffs == n_mm, f"Expected {n_mm} mismatches, got {diffs}"

    def test_skips_impossible_mismatch(self):
        short = "AC"
        results = generate_mismatch_sequences(short, mismatch_counts=(5,), num_per_count=1)
        assert len(results) == 0


class TestGenerateOfftargetSequences:
    """Test legacy generate_offtarget_sequences."""

    def test_returns_list(self):
        offtargets = generate_offtarget_sequences(TARGET_REGION, num_scrambled=1, num_mismatch=1, seed=1)
        assert isinstance(offtargets, list)
        assert len(offtargets) >= 1

    def test_structure_preserved(self):
        offtargets = generate_offtarget_sequences(TARGET_REGION, num_scrambled=1, num_mismatch=1, seed=1)
        for ot in offtargets:
            assert ot.startswith("AAAAAA") and ot.endswith("AAAAAA")
            assert len(ot) == len(DUMMY_TARGET_RNA)


class TestGenerateEvaluationJsons:
    """Test generate_evaluation_jsons. No subprocess."""

    def test_creates_off_and_on_jsons(self, variant_fasta, sample_metadata_json, tmpdir):
        off_path, on_path = generate_evaluation_jsons(
            variant_fasta, "test_cas13", sample_metadata_json, str(tmpdir)
        )
        assert Path(off_path).exists()
        assert Path(on_path).exists()
        off = json.loads(Path(off_path).read_text())
        on = json.loads(Path(on_path).read_text())
        assert len(off[0]["sequences"]) == 2  # protein + crRNA
        assert len(on[0]["sequences"]) == 3   # protein + crRNA + target
        assert "protein" in str(off)
        assert "rna" in str(off)

    def test_crrna_lookup_override(self, variant_fasta, sample_metadata_json, tmpdir):
        off_path, on_path = generate_evaluation_jsons(
            variant_fasta, "other_id", sample_metadata_json, str(tmpdir), crrna_lookup_id="test_cas13"
        )
        assert Path(off_path).exists()

    def test_missing_metadata_raises(self, variant_fasta, tmpdir):
        with pytest.raises(ValueError, match="not found"):
            generate_evaluation_jsons(variant_fasta, "nonexistent_id", "/no/file.json", str(tmpdir))


class TestGenerateOfftargetJson:
    """Test generate_offtarget_json."""

    def test_creates_valid_json(self, variant_fasta, sample_metadata_json, tmpdir):
        off_rna = "AAAAAA" + "C" * 24 + "AAAAAA"  # all mismatches
        path = generate_offtarget_json(
            variant_fasta, "test_cas13", sample_metadata_json, off_rna, str(tmpdir), suffix="0"
        )
        assert Path(path).exists()
        data = json.loads(Path(path).read_text())
        assert len(data[0]["sequences"]) == 3
        assert "offtarget" in data[0]["name"]

    def test_missing_baseline_raises(self, variant_fasta, tmpdir):
        meta = tmpdir / "meta.json"
        meta.write_text("{}")
        with pytest.raises(ValueError, match="not found"):
            generate_offtarget_json(variant_fasta, "x", str(meta), "A" * 36, str(tmpdir))
