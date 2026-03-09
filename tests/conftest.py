"""
CASCADE test fixtures. All tests run without GPU/Protenix/PXDesign.
Run from project root: pytest tests/ -v
"""
import json
import os
import sys
from pathlib import Path

# Add scripts/ to path for imports
PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import tempfile

import pytest


# --- Minimal PDB with two residues for HEPN distance testing ---
# Residues at positions 10 and 50, CA atoms at (0,0,0) and (30,0,0) = 30 Angstroms
MINIMAL_PDB = """ATOM      1  N   ALA A  10       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A  10       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A  10       0.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   ALA A  10       0.000   0.000   0.000  1.00  0.00           O
ATOM      5  N   ALA A  11       0.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  ALA A  11       0.000   0.000   0.000  1.00  0.00           C
ATOM      7  N   HIS A  50      30.000   0.000   0.000  1.00  0.00           N
ATOM      8  CA  HIS A  50      30.000   0.000   0.000  1.00  0.00           C
ATOM      9  C   HIS A  50      30.000   0.000   0.000  1.00  0.00           C
ATOM     10  O   HIS A  50      30.000   0.000   0.000  1.00  0.00           O
END
"""

# PDB with 15 Angstrom distance
MINIMAL_PDB_15A = """ATOM      1  CA  ALA A  10       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  HIS A  50      15.000   0.000   0.000  1.00  0.00           C
END
"""

# --- Biologically realistic test sequences ---
# These sequences must have TWO R.{3,6}H motifs separated by ≥150 residues
# (matching the MIN_HEPN_SEPARATION filter in identify_hepn_domains and
# get_catalytic_histidine_indices).
#
# Layout (0-based): HEPN1 motif at ~pos 6 (RAILXH), HEPN2 motif at ~pos 200 (RVVVXH)
# Separation: 200 - 6 = 194 residues → passes MIN_HEPN_SEPARATION (150)
# Total length: ~250 residues

_PRE_HEPN1 = "MAAAAA"           # 6 residues (NTD region)
_HEPN1_MOTIF = "RAILXH"         # 6 residues (HEPN1 catalytic motif, R at pos 6)
_INTER_DOMAIN = "G" * 188       # 188 residues (inter-domain: Helical-2 + IDL + HEPN2-upstream)
_HEPN2_MOTIF = "RVVVXH"         # 6 residues (HEPN2 catalytic motif, R at pos 200)
_C_TERMINAL = "G" * 44          # 44 residues (C-terminal tail)

_TEST_SEQ = _PRE_HEPN1 + _HEPN1_MOTIF + _INTER_DOMAIN + _HEPN2_MOTIF + _C_TERMINAL
assert len(_TEST_SEQ) == 250, f"Expected 250, got {len(_TEST_SEQ)}"

# Variant with G->P at position 13 (1-based), which is in the inter-domain region
_VARIANT_SEQ = _TEST_SEQ[:12] + "P" + _TEST_SEQ[13:]

FASTA_WITH_HEPN = f">test_cas13\n{_TEST_SEQ}\n"
FASTA_VARIANT_SINGLE_MUT = f">test_cas13_variant\n{_VARIANT_SEQ}\n"
FASTA_BASELINE = f">baseline_id\n{_TEST_SEQ}\n"

# Sequence with only one HEPN motif (should fail)
FASTA_ONE_HEPN = ">bad_seq\nMAAAAARAILXHGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"

# Sequence with no HEPN motif
FASTA_NO_HEPN = ">no_hepn\nMGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"

# Sequence with two motifs but too close together (< 150 residues apart)
FASTA_CLOSE_HEPN = ">close_hepn\nMAAAAARAILXH" + "G" * 100 + "RVVVXH" + "G" * 40 + "\n"


@pytest.fixture
def tmpdir():
    """Temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


@pytest.fixture
def data_dir(tmpdir):
    """data/mined_hits style directory with FASTA and CSV."""
    d = tmpdir / "data" / "mined_hits"
    d.mkdir(parents=True)
    return d


@pytest.fixture
def sample_fasta(data_dir):
    """Sample FASTA with HEPN motifs (250 residues, HEPN1 at 6, HEPN2 at 200)."""
    p = data_dir / "test_hits.fasta"
    p.write_text(FASTA_WITH_HEPN)
    return str(p)


@pytest.fixture
def sample_csv(data_dir):
    """Sample metadata CSV matching FASTA IDs."""
    p = data_dir / "test_hits_metadata.csv"
    content = """sequence_id,repeat_domains,sra_accession,score
test_cas13,AAACCCGGGTTT,SRR123,0.9
"""
    p.write_text(content)
    return str(p)


@pytest.fixture
def sample_metadata_json(tmpdir):
    """Variant domain metadata as produced by 01_parse_and_annotate.
    Uses expanded HEPN boundaries: -60 upstream, +100 downstream of catalytic R."""
    p = tmpdir / "variant_domain_metadata.json"
    # HEPN1 R is at 0-based pos 6; domain = [max(0,6-60), 6+100] = [0, 106]
    # HEPN2 R is at 0-based pos 200; domain = [max(106,200-60), 200+100] = [140, 300]
    meta = {
        "test_cas13": {
            "sequence_length": 250,
            "domains": {
                "HEPN1": {"start": 0, "end": 106},
                "HEPN2": {"start": 140, "end": 300},
            },
            "crRNA_repeat_used": "AAACCCGGGUUU",
        },
        "baseline_id": {
            "sequence_length": 250,
            "domains": {
                "HEPN1": {"start": 0, "end": 106},
                "HEPN2": {"start": 140, "end": 300},
            },
            "crRNA_repeat_used": "AAACCCGGGUUU",
        },
    }
    p.write_text(json.dumps(meta, indent=2))
    return str(p)


@pytest.fixture
def minimal_pdb(tmpdir):
    """Minimal PDB with two CA atoms 30A apart at residues 10 and 50."""
    p = tmpdir / "model.pdb"
    p.write_text(MINIMAL_PDB)
    return str(p)


@pytest.fixture
def minimal_pdb_15a(tmpdir):
    """Minimal PDB with 15A distance."""
    p = tmpdir / "model_15a.pdb"
    p.write_text(MINIMAL_PDB_15A)
    return str(p)


@pytest.fixture
def protenix_summary_json(tmpdir):
    """Mock Protenix output_summary.json."""
    p = tmpdir / "model_summary.json"
    data = {
        "iptm": 0.92,
        "ptm": 0.88,
        "ranking_score": 0.9,
        "af2_ig": 0.85,
    }
    p.write_text(json.dumps(data))
    return str(p)


@pytest.fixture
def variant_fasta(tmpdir):
    """Variant FASTA for evaluation JSON generation."""
    p = tmpdir / "variant_001.fasta"
    p.write_text(FASTA_VARIANT_SINGLE_MUT)
    return str(p)


@pytest.fixture
def baseline_json(tmpdir, sample_metadata_json):
    """Baseline Protenix JSON (from jsons/)."""
    meta = json.loads(Path(sample_metadata_json).read_text())
    crrna = meta["baseline_id"]["crRNA_repeat_used"] + "GUCGACUGACGUACGUACGUACGU"
    # Protenix format: proteinChain/rnaSequence
    payload = [{
        "name": "baseline_id",
        "sequences": [
            {"proteinChain": {"sequence": _TEST_SEQ, "count": 1}},
            {"rnaSequence": {"sequence": crrna, "count": 1}},
            {"rnaSequence": {"sequence": "AAAAAAACGUACGUACGUACGUCAGUCGACAAAAAA", "count": 1}},
        ],
    }]
    d = tmpdir / "jsons"
    d.mkdir()
    p = d / "baseline_id.json"
    p.write_text(json.dumps(payload, indent=2))
    return str(p)
