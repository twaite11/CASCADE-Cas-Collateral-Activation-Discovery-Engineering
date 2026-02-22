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

# Sequence with R.{4,6}H motif: RAILXH and RVVVXH (1-based His at 12 and 37)
FASTA_WITH_HEPN = """>test_cas13
MAAAAARAILXHGGGGGGGGGGGGGGGGGGGRVVVXHGGGGGGGG
"""

# Variant with G->P at position 13 (1-based)
FASTA_VARIANT_SINGLE_MUT = """>test_cas13_variant
MAAAAARAILXHPGGGGGGGGGGGGGGGGGGRVVVXHGGGGGGGG
"""

# Baseline for mutation comparison
FASTA_BASELINE = """>baseline_id
MAAAAARAILXHGGGGGGGGGGGGGGGGGGGRVVVXHGGGGGGGG
"""

# Sequence with only one HEPN motif (should fail)
FASTA_ONE_HEPN = """>bad_seq
MAAAAARAILXHGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
"""

# Sequence with no HEPN motif
FASTA_NO_HEPN = """>no_hepn
MGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
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
    """Variant domain metadata as produced by 01_parse_and_annotate."""
    p = tmpdir / "variant_domain_metadata.json"
    meta = {
        "test_cas13": {
            "sequence_length": 45,
            "domains": {
                "HEPN1": {"start": 0, "end": 85},
                "HEPN2": {"start": 20, "end": 105},
            },
            "crRNA_repeat_used": "AAACCCGGGUUU",
        },
        "baseline_id": {
            "sequence_length": 45,
            "domains": {
                "HEPN1": {"start": 0, "end": 85},
                "HEPN2": {"start": 20, "end": 105},
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
            {"proteinChain": {"sequence": "MAAAAARAILXHGGGGGGGGGGGGGGGGGGGRVVVXHGGGGGGGG", "count": 1}},
            {"rnaSequence": {"sequence": crrna, "count": 1}},
            {"rnaSequence": {"sequence": "AAAAAAACGUACGUACGUACGUCAGUCGACAAAAAA", "count": 1}},
        ],
    }]
    d = tmpdir / "jsons"
    d.mkdir()
    p = d / "baseline_id.json"
    p.write_text(json.dumps(payload, indent=2))
    return str(p)
