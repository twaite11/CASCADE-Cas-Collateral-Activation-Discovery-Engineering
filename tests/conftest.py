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


# --- Minimal PDB with two HIS residues for HEPN distance testing ---
# Residue 10 (HIS): CA at (0,0,0), NE2 at (2,0,0)
# Residue 50 (HIS): CA at (30,0,0), NE2 at (28,0,0)
# CA-CA = 30 A, NE2-NE2 = 26 A
MINIMAL_PDB = """ATOM      1  N   HIS A  10       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  HIS A  10       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   HIS A  10       0.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   HIS A  10       0.000   0.000   0.000  1.00  0.00           O
ATOM      5  CG  HIS A  10       1.000   0.000   0.000  1.00  0.00           C
ATOM      6  ND1 HIS A  10       1.500   0.000   0.000  1.00  0.00           N
ATOM      7  NE2 HIS A  10       2.000   0.000   0.000  1.00  0.00           N
ATOM      8  N   ALA A  11       0.000   0.000   0.000  1.00  0.00           N
ATOM      9  CA  ALA A  11       0.000   0.000   0.000  1.00  0.00           C
ATOM     10  N   HIS A  50      30.000   0.000   0.000  1.00  0.00           N
ATOM     11  CA  HIS A  50      30.000   0.000   0.000  1.00  0.00           C
ATOM     12  C   HIS A  50      30.000   0.000   0.000  1.00  0.00           C
ATOM     13  O   HIS A  50      30.000   0.000   0.000  1.00  0.00           O
ATOM     14  CG  HIS A  50      29.000   0.000   0.000  1.00  0.00           C
ATOM     15  ND1 HIS A  50      28.500   0.000   0.000  1.00  0.00           N
ATOM     16  NE2 HIS A  50      28.000   0.000   0.000  1.00  0.00           N
END
"""

# PDB with NE2-NE2 = 11 A (CA-CA = 15 A)
# Residue 10 (HIS): CA at (0,0,0), NE2 at (2,0,0)
# Residue 50 (HIS): CA at (15,0,0), NE2 at (13,0,0)
MINIMAL_PDB_15A = """ATOM      1  CA  HIS A  10       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  NE2 HIS A  10       2.000   0.000   0.000  1.00  0.00           N
ATOM      3  CA  HIS A  50      15.000   0.000   0.000  1.00  0.00           C
ATOM      4  NE2 HIS A  50      13.000   0.000   0.000  1.00  0.00           N
END
"""

# PDB with only CA atoms (no side-chain) — tests fallback to CA
MINIMAL_PDB_CA_ONLY = """ATOM      1  N   HIS A  10       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  HIS A  10       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  N   HIS A  50      30.000   0.000   0.000  1.00  0.00           N
ATOM      4  CA  HIS A  50      30.000   0.000   0.000  1.00  0.00           C
END
"""

# Sequence with two R.{4,6}H HEPN motifs separated by ~300 residues.
# HEPN1 motif at pos ~110, HEPN2 motif at ~410.  Total length = 500.
_PAD_100 = "G" * 100
_PAD_200 = "G" * 200
_HEPN_SEQ = f"M{'A' * 109}RAILXH{_PAD_200}{'A' * 94}RVVVXH{_PAD_100}"  # 500 aa
FASTA_WITH_HEPN = f""">test_cas13
{_HEPN_SEQ}
"""

# Variant with A->P at position 113 (1-based, just after first HEPN motif)
_VARIANT_SEQ = _HEPN_SEQ[:112] + "P" + _HEPN_SEQ[113:]
FASTA_VARIANT_SINGLE_MUT = f""">test_cas13_variant
{_VARIANT_SEQ}
"""

FASTA_BASELINE = f""">baseline_id
{_HEPN_SEQ}
"""

# Sequence with only one HEPN motif (should fail)
FASTA_ONE_HEPN = f""">bad_seq
M{'A' * 109}RAILXH{'G' * 390}
"""

# Sequence with no HEPN motif
FASTA_NO_HEPN = f""">no_hepn
M{'G' * 499}
"""


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
    """Sample FASTA with HEPN motifs."""
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
    HEPN1 motif at ~110, HEPN2 motif at ~410 in 500-residue sequences."""
    p = tmpdir / "variant_domain_metadata.json"
    meta = {
        "test_cas13": {
            "sequence_length": 500,
            "domains": {
                "HEPN1": {"start": 80, "end": 190},
                "HEPN2": {"start": 380, "end": 490},
            },
            "crRNA_repeat_used": "AAACCCGGGUUU",
        },
        "baseline_id": {
            "sequence_length": 500,
            "domains": {
                "HEPN1": {"start": 80, "end": 190},
                "HEPN2": {"start": 380, "end": 490},
            },
            "crRNA_repeat_used": "AAACCCGGGUUU",
        },
    }
    p.write_text(json.dumps(meta, indent=2))
    return str(p)


@pytest.fixture
def minimal_pdb(tmpdir):
    """Minimal PDB with two HIS residues: NE2-NE2 = 26 A, CA-CA = 30 A."""
    p = tmpdir / "model.pdb"
    p.write_text(MINIMAL_PDB)
    return str(p)


@pytest.fixture
def minimal_pdb_15a(tmpdir):
    """Minimal PDB with two HIS residues: NE2-NE2 = 11 A, CA-CA = 15 A."""
    p = tmpdir / "model_15a.pdb"
    p.write_text(MINIMAL_PDB_15A)
    return str(p)


@pytest.fixture
def minimal_pdb_ca_only(tmpdir):
    """Minimal PDB with HIS residues but only CA atoms (no side-chain)."""
    p = tmpdir / "model_ca_only.pdb"
    p.write_text(MINIMAL_PDB_CA_ONLY)
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
    payload = [{
        "name": "baseline_id",
        "sequences": [
            {"proteinChain": {"sequence": _HEPN_SEQ, "count": 1}},
            {"rnaSequence": {"sequence": crrna, "count": 1}},
            {"rnaSequence": {"sequence": "AAAAAAACGUACGUACGUACGUCAGUCGACAAAAAA", "count": 1}},
        ],
    }]
    d = tmpdir / "jsons"
    d.mkdir()
    p = d / "baseline_id.json"
    p.write_text(json.dumps(payload, indent=2))
    return str(p)
