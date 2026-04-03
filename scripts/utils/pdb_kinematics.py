import json
import os
import glob
import shutil
import subprocess
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO
import warnings
from Bio import BiopythonWarning

# Suppress minor PDB format warnings for cleaner RunPod logs
warnings.simplefilter('ignore', BiopythonWarning)

_RUST_STRUCTSCORE_BIN = shutil.which("cascade_structscore")
if not _RUST_STRUCTSCORE_BIN:
    _candidate = os.path.join(os.path.dirname(__file__), "..", "..", "rust", "target", "release", "cascade_structscore")
    if os.name == "nt":
        _candidate += ".exe"
    if os.path.isfile(_candidate):
        _RUST_STRUCTSCORE_BIN = _candidate


def _run_rust_structscore(args, timeout=60):
    """Run cascade_structscore with args. Returns stdout string on success, None on failure."""
    if not _RUST_STRUCTSCORE_BIN:
        return None
    try:
        r = subprocess.run(
            [_RUST_STRUCTSCORE_BIN] + args,
            capture_output=True, text=True, timeout=timeout,
        )
        if r.returncode == 0:
            return r.stdout.strip()
    except Exception:
        pass
    return None


def find_structure_files(base_dir):
    """
    Search recursively for structure files, preferring CIF over PDB.
    Returns list of paths (CIF first if any, else PDB).
    """
    if not base_dir or not os.path.isdir(base_dir):
        return []
    rust_out = _run_rust_structscore(["find-structures", "--dir", str(base_dir)])
    if rust_out is not None:
        try:
            return json.loads(rust_out)
        except (json.JSONDecodeError, TypeError):
            pass
    cifs = sorted(glob.glob(os.path.join(base_dir, "**", "*.cif"), recursive=True))
    if cifs:
        return cifs
    pdbs = sorted(glob.glob(os.path.join(base_dir, "**", "*.pdb"), recursive=True))
    return pdbs


def _load_structure(path):
    """Load a structure from PDB or CIF; returns Biopython Structure object."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("structure", path)

def extract_protenix_scores(summary_json_path):
    """
    Parses the Protenix _summary.json output to extract the true ipTM, pTM, and AF2-IG scores.
    """
    if not os.path.exists(summary_json_path):
        raise FileNotFoundError(f"Protenix summary file not found: {summary_json_path}")

    rust_out = _run_rust_structscore(["extract-scores", "--summary", str(summary_json_path)])
    if rust_out is not None:
        try:
            return json.loads(rust_out)
        except (json.JSONDecodeError, TypeError):
            pass

    with open(summary_json_path, 'r') as f:
        data = json.load(f)
        
    # Standard AlphaFold3/Protenix summary dictionary keys
    iptm = float(data.get('iptm', 0.0))
    ptm = float(data.get('ptm', 0.0))
    ranking_score = float(data.get('ranking_score', 0.0))
    
    # Extract AF2-IG (Interface Gap/Confidence) score which is critical for scoring RNP complexes
    af2_ig = float(data.get('af2_ig', data.get('af2_ig_score', 0.0)))
    
    return {
        "iptm": iptm,
        "ptm": ptm,
        "ranking_score": ranking_score,
        "af2_ig": af2_ig
    }

def calculate_hepn_shift(structure_path, hepn1_his_idx, hepn2_his_idx, protein_chain_id='A'):
    """
    Parses a PDB or CIF file and calculates the 3D Euclidean distance (in Angstroms)
    between the Alpha-Carbons of the two catalytic Histidines in the HEPN domains.
    Accepts .pdb or .cif; uses Biopython MMCIFParser for CIF.
    """
    if not os.path.exists(structure_path):
        raise FileNotFoundError(f"Structure file not found: {structure_path}")

    rust_out = _run_rust_structscore([
        "hepn-distance",
        "--structure", str(structure_path),
        "--h1-idx", str(hepn1_his_idx),
        "--h2-idx", str(hepn2_his_idx),
        "--chain", str(protein_chain_id),
    ])
    if rust_out is not None:
        try:
            return float(rust_out)
        except (ValueError, TypeError):
            pass

    structure = _load_structure(structure_path)
    
    try:
        model = structure[0]
        chain = model[protein_chain_id]
        
        res1 = chain[hepn1_his_idx]
        res2 = chain[hepn2_his_idx]
        
        coord1 = res1['CA'].get_coord()
        coord2 = res2['CA'].get_coord()
        
        distance_angstroms = np.linalg.norm(coord1 - coord2)
        return float(distance_angstroms)
        
    except KeyError as e:
        raise ValueError(f"Could not find required residue or chain in structure {structure_path}. Error: {e}")


def cif_to_pdb(cif_path, pdb_path=None):
    """
    Convert a CIF file to PDB using Biopython.
    Returns path to the written PDB file.
    """
    if not os.path.exists(cif_path):
        raise FileNotFoundError(f"CIF file not found: {cif_path}")
    if pdb_path is None:
        pdb_path = cif_path.replace(".cif", ".pdb").replace(".CIF", ".pdb")
    structure = _load_structure(cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)
    return pdb_path