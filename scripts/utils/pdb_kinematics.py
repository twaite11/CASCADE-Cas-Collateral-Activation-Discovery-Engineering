import json
import logging
import os
import glob
import shutil
import subprocess
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO
import warnings
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

log = logging.getLogger(__name__)

_HIS_RESNAMES = frozenset(("HIS", "HSE", "HSD", "HSP"))

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

def _compute_interface_score_from_pae(data):
    """
    Derive an AF2-style interface/interaction-geometry score from Protenix
    chain-pair metrics when the explicit af2_ig key is absent.

    Strategy (in priority order):
    1. chain_pair_iptm matrix: mean of off-diagonal elements gives cross-chain
       interface confidence (protein-RNA interaction quality).
    2. Weighted combination: 0.8*iptm + 0.2*ptm as a proxy when no
       chain_pair data is available but global scores exist.
    """
    cp_iptm = data.get("chain_pair_iptm")
    if cp_iptm and isinstance(cp_iptm, (list, tuple)):
        off_diag = []
        n = len(cp_iptm)
        for i in range(n):
            row = cp_iptm[i]
            if not isinstance(row, (list, tuple)):
                continue
            for j in range(len(row)):
                if i != j:
                    try:
                        off_diag.append(float(row[j]))
                    except (ValueError, TypeError):
                        pass
        if off_diag:
            return float(np.mean(off_diag))

    iptm = float(data.get("iptm", 0.0) or 0.0)
    ptm = float(data.get("ptm", 0.0) or 0.0)
    if iptm > 0.0 or ptm > 0.0:
        return 0.8 * iptm + 0.2 * ptm

    return 0.0


def extract_protenix_scores(summary_json_path):
    """
    Parses the confidence/summary JSON output to extract ipTM, pTM, and AF2-IG scores.
    Compatible with both Protenix (*_summary*.json) and cattle-prod (*_summary_confidence.json).

    AF2-IG derivation: if the explicit af2_ig/af2_ig_score key is absent or zero,
    computes an interface score from chain_pair_iptm (off-diagonal mean) or falls
    back to a weighted iptm/ptm proxy.
    """
    if not os.path.exists(summary_json_path):
        raise FileNotFoundError(f"Summary/confidence file not found: {summary_json_path}")

    rust_out = _run_rust_structscore(["extract-scores", "--summary", str(summary_json_path)])
    if rust_out is not None:
        try:
            parsed = json.loads(rust_out)
            if float(parsed.get("af2_ig", 0.0)) > 0.0:
                return parsed
        except (json.JSONDecodeError, TypeError, ValueError):
            pass

    with open(summary_json_path, 'r') as f:
        data = json.load(f)
        
    iptm = float(data.get('iptm', 0.0) or 0.0)
    ptm = float(data.get('ptm', 0.0) or 0.0)
    ranking_score = float(data.get('ranking_score', 0.0) or 0.0)
    af2_ig = float(data.get('af2_ig', data.get('af2_ig_score', 0.0)) or 0.0)

    if af2_ig == 0.0:
        af2_ig = _compute_interface_score_from_pae(data)
    
    return {
        "iptm": iptm,
        "ptm": ptm,
        "ranking_score": ranking_score,
        "af2_ig": af2_ig
    }

def _get_catalytic_coord(residue, label=""):
    """
    Get the most catalytically-relevant coordinate for a HEPN residue.

    For histidines the NE2 atom is the functional group that coordinates Mg2+
    and the RNA substrate. Falls back through ND1 -> CG -> CA if side-chain
    atoms are missing (e.g. backbone-only predictions).
    """
    resname = residue.get_resname().strip().upper()
    if resname not in _HIS_RESNAMES:
        log.warning(
            f"{label} catalytic residue is {resname} (id={residue.id}), not HIS — "
            f"HEPN motif may have anchored incorrectly"
        )

    for atom_name in ('NE2', 'ND1', 'CG', 'CA'):
        if atom_name in residue:
            return residue[atom_name].get_coord()
    raise KeyError(f"No suitable atom found in {resname} {residue.id}")


def calculate_hepn_shift(structure_path, hepn1_his_idx, hepn2_his_idx, protein_chain_id='A'):
    """
    Parses a PDB or CIF file and calculates the 3D Euclidean distance (in Angstroms)
    between the catalytic histidines in the two HEPN domains.

    Prefers the side-chain NE2 atom (the functional imidazole nitrogen that
    coordinates Mg2+ and the RNA phosphodiester backbone). Falls back to
    ND1 -> CG -> CA when side-chain atoms are absent.
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
        
        std_residues = [r for r in chain.get_residues() if r.id[0] == ' ']
        idx1 = hepn1_his_idx - 1
        idx2 = hepn2_his_idx - 1
        if 0 <= idx1 < len(std_residues) and 0 <= idx2 < len(std_residues):
            res1, res2 = std_residues[idx1], std_residues[idx2]
        else:
            res1 = chain[hepn1_his_idx]
            res2 = chain[hepn2_his_idx]
        
        coord1 = _get_catalytic_coord(res1, label="HEPN1")
        coord2 = _get_catalytic_coord(res2, label="HEPN2")
        
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