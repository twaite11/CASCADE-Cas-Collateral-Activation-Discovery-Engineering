import json
import os
import numpy as np
from Bio.PDB import PDBParser
import warnings
from Bio import BiopythonWarning

# Suppress minor PDB format warnings for cleaner RunPod logs
warnings.simplefilter('ignore', BiopythonWarning)

def extract_protenix_scores(summary_json_path):
    """
    Parses the Protenix _summary.json output to extract the true ipTM, pTM, and AF2-IG scores.
    """
    if not os.path.exists(summary_json_path):
        raise FileNotFoundError(f"Protenix summary file not found: {summary_json_path}")
        
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

def calculate_hepn_shift(pdb_path, hepn1_his_idx, hepn2_his_idx, protein_chain_id='A'):
    """
    Parses a PDB file and calculates the 3D Euclidean distance (in Angstroms) 
    between the Alpha-Carbons of the two catalytic Histidines in the HEPN domains.
    """
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("cas13_complex", pdb_path)
    
    try:
        # Navigate the PDB hierarchy: Structure -> Model (0) -> Chain -> Residue -> Atom
        model = structure[0]
        chain = model[protein_chain_id]
        
        # Biopython parses 1-based PDB residue indices. Ensure your SQLite indices match.
        res1 = chain[hepn1_his_idx]
        res2 = chain[hepn2_his_idx]
        
        # Extract the 3D numpy coordinate arrays for the Alpha-Carbons (CA)
        coord1 = res1['CA'].get_coord()
        coord2 = res2['CA'].get_coord()
        
        # Calculate the Euclidean distance
        distance_angstroms = np.linalg.norm(coord1 - coord2)
        return float(distance_angstroms)
        
    except KeyError as e:
        raise ValueError(f"Could not find required residue or chain in PDB {pdb_path}. Error: {e}")