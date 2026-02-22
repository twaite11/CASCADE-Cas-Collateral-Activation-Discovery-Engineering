"""
PXDesign wrapper for CASCADE evolution.
Uses the actual PXDesign CLI: pxdesign infer -i <yaml> -o <dir> --N_sample N
Generation only; we compute Protenix scores downstream. Generates YAML from baseline structure + metadata,
runs inference, parses CIF output to variant FASTAs.
"""
import subprocess
import json
import os
import glob
import logging
from pathlib import Path

log = logging.getLogger(__name__)

# Map 3-letter to 1-letter amino acid codes (standard + common)
_AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _get_structure_chain_ids(structure_path: str):
    """Return list of chain IDs in the structure (e.g. ['A','B','C'] for ternary)."""
    from Bio.PDB import MMCIFParser, PDBParser
    ext = os.path.splitext(structure_path)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext == ".cif" else PDBParser(QUIET=True)
    struct = parser.get_structure("s", structure_path)
    return [c.id for c in struct[0].get_chains()]


def _sequence_from_structure(structure_path: str, chain_id: str = "A") -> str:
    """Extract protein sequence from CIF or PDB using Biopython."""
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.PDB.Polypeptide import three_to_one
    ext = os.path.splitext(structure_path)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext == ".cif" else PDBParser(QUIET=True)
    struct = parser.get_structure("s", structure_path)
    chain = struct[0][chain_id]
    seq = []
    for res in chain:
        if res.id[0] != " ":
            continue
        resname = res.get_resname()
        try:
            seq.append(three_to_one(resname))
        except KeyError:
            seq.append(_AA3_TO_1.get(resname, "X"))
    return "".join(seq)


def _sequence_from_structure_last_chain(structure_path: str) -> str:
    """Extract sequence from the last chain (PXDesign outputs binder as final chain)."""
    chain_ids = _get_structure_chain_ids(structure_path)
    if not chain_ids:
        return ""
    return _sequence_from_structure(structure_path, chain_ids[-1])


def _sequence_from_fasta_or_json(fasta_path: str, json_path: str, variant_id: str) -> str:
    """Get baseline sequence from FASTA or base JSON."""
    if fasta_path and os.path.exists(fasta_path):
        with open(fasta_path) as f:
            return "".join(l.strip() for l in f if not l.startswith(">"))
    if os.path.exists(json_path):
        with open(json_path) as f:
            data = json.load(f)
        ent = data[0]["sequences"][0]
        prot = ent.get("proteinChain", ent.get("protein", {}))
        return prot.get("sequence", "")
    return ""


def generate_frozen_rec_config(metadata_db_path, variant_id, metadata_override=None):
    """
    Extracts boundaries for REC, HEPN1, and HEPN2 domains.
    Returns a coords dict with rec_end, hepn1/hepn2 boundaries, seq_len, binder_length (linkers only).
    binder_length = linker1 + linker2 so PXDesign designs a smaller, tighter mass.
    """
    if metadata_override and variant_id in metadata_override:
        variant_data = metadata_override[variant_id]
    else:
        with open(metadata_db_path, 'r') as f:
            metadata = json.load(f)
        variant_data = metadata.get(variant_id)
        if not variant_data:
            raise ValueError(f"Variant {variant_id} not found in metadata.")

    hepn1 = variant_data["domains"]["HEPN1"]
    hepn2 = variant_data["domains"]["HEPN2"]
    hepn1_start, hepn1_end = hepn1["start"], hepn1["end"]
    hepn2_start, hepn2_end = hepn2["start"], hepn2["end"]
    rec_end = max(1, hepn1_start - 10)
    seq_len = variant_data.get("sequence_length", hepn2_end + 20)

    # binder_length = linkers only (REC→HEPN1 and HEPN1→HEPN2); avoids 600-aa de novo design
    linker1_len = max(0, hepn1_start - rec_end - 1)
    linker2_len = max(0, hepn2_start - hepn1_end - 1)
    binder_length = linker1_len + linker2_len
    MIN_BINDER_LENGTH = 20  # PXDesign may struggle with very short binders
    if binder_length < MIN_BINDER_LENGTH:
        binder_length = MIN_BINDER_LENGTH

    return {
        "rec_end": rec_end,
        "hepn1_start": hepn1_start, "hepn1_end": hepn1_end,
        "hepn2_start": hepn2_start, "hepn2_end": hepn2_end,
        "seq_len": seq_len,
        "binder_length": binder_length,
        "linker1_len": linker1_len,
        "linker2_len": linker2_len,
    }


def _build_pxdesign_yaml(
    baseline_structure: str,
    variant_id: str,
    coords: dict,
    output_dir: str,
) -> str:
    """Build PXDesign YAML: REC crop, RNA chains B/C as fixed context, linker-only binder_length."""
    yaml_path = os.path.join(output_dir, f"{variant_id}_pxdesign_input.yaml")
    abs_structure = os.path.abspath(baseline_structure)
    if not os.path.exists(abs_structure):
        raise FileNotFoundError(f"Baseline structure not found: {abs_structure}")

    rec_end = coords["rec_end"]
    binder_length = coords["binder_length"]
    hotspots = list(range(max(1, rec_end - 5), rec_end + 1))
    if not hotspots:
        hotspots = [rec_end]

    chains = {
        "A": {"crop": [f"1-{rec_end}"], "hotspots": hotspots[:10]},
    }
    chain_ids = _get_structure_chain_ids(abs_structure)
    if "B" in chain_ids:
        chains["B"] = "all"
    if "C" in chain_ids:
        chains["C"] = "all"

    cfg = {
        "task_name": variant_id,
        "binder_length": binder_length,
        "target": {"file": abs_structure, "chains": chains},
    }
    try:
        import yaml
        with open(yaml_path, "w") as f:
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
    except ImportError:
        import json
        with open(yaml_path.replace(".yaml", ".json"), "w") as f:
            json.dump(cfg, f, indent=2)
        yaml_path = yaml_path.replace(".yaml", ".json")
    return yaml_path


def run_pxdesign_generation(
    baseline_structure: str,
    variant_id: str,
    metadata_path: str,
    bias_json_path: str,
    output_dir: str,
    variant_count: int = 50,
    metadata_override=None,
    baseline_fasta_path: str = None,
    base_json_dir: str = None,
):
    """
    Runs PXDesign via: pxdesign infer -i <yaml> -o <dir> --N_sample N
    Generation only (no evaluation); we compute scores ourselves with Protenix later.
    Designs linkers only; stitches wild-type HEPN1/HEPN2 into output.
    Fallback: returns baseline as variant_fallback_* when stitching fails (orchestrator applies penalty).
    """
    from utils.hepn_structural_stitch import stitch_hepn_into_binder

    log.info(f"[PXDesign] Generating {variant_count} designs for {variant_id} (5-30 min)...")
    os.makedirs(output_dir, exist_ok=True)

    coords = generate_frozen_rec_config(metadata_path, variant_id, metadata_override)
    rec_end = coords["rec_end"]
    binder_length = coords["binder_length"]

    if base_json_dir:
        base_json = os.path.join(base_json_dir, f"{variant_id}.json")
    else:
        meta_dir = os.path.dirname(os.path.abspath(metadata_path))
        base_json = os.path.join(meta_dir, "..", "jsons", f"{variant_id}.json")
    full_seq = _sequence_from_fasta_or_json(baseline_fasta_path, base_json, variant_id)
    if not full_seq:
        full_seq = _sequence_from_structure(baseline_structure)
    if len(full_seq) < coords["hepn2_end"]:
        log.warning(f"Baseline sequence shorter than HEPN2 end; stitching may fail")

    yaml_path = _build_pxdesign_yaml(baseline_structure, variant_id, coords, output_dir)
    abs_out = os.path.abspath(output_dir)

    cmd = [
        "pxdesign", "infer",
        "-i", yaml_path,
        "-o", abs_out,
        "--N_sample", str(variant_count),
        "--dtype", "bf16",
    ]
    if bias_json_path and os.path.exists(bias_json_path):
        log.info(f"  -> RL bias matrix available but PXDesign infer does not support --bias_aa_json; skipping")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            timeout=3600,
            cwd=os.path.dirname(os.path.abspath(__file__)) or ".",
        )
    except subprocess.CalledProcessError as e:
        log.error(f"PXDesign failed:\n{e.stderr}")
        raise
    except subprocess.TimeoutExpired:
        log.error("PXDesign timed out after 2 hours")
        raise

    def _write_variant_or_fallback(i: int, binder_seq: str) -> str:
        """Stitch HEPN into binder; on failure, write baseline as fallback (orchestrator applies penalty)."""
        full = stitch_hepn_into_binder(binder_seq, full_seq, coords)
        if full:
            name = f"{variant_id}_variant_{i}"
            fasta_path = os.path.join(output_dir, f"{name}.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">{name}\n{full}\n")
            return fasta_path
        # Fallback: baseline sequence; name contains "fallback" so orchestrator applies penalty
        fallback_path = os.path.join(output_dir, f"{variant_id}_variant_fallback_{i}.fasta")
        with open(fallback_path, "w") as f:
            f.write(f">{variant_id}_variant_fallback_{i}\n{full_seq}\n")
        log.warning(f"Stitching failed for design {i}; writing baseline as fallback (will receive penalty)")
        return fallback_path

    design_out = os.path.join(abs_out, "design_outputs")
    summary_csvs = glob.glob(os.path.join(design_out, "**", "summary.csv"), recursive=True)
    if not summary_csvs:
        log.warning("No summary.csv found; checking orig_designed and predictions")
        pred_glob = glob.glob(os.path.join(abs_out, "**", "predictions", "*.cif"), recursive=True)
        if not pred_glob:
            return []
        fastas = []
        for i, cif_path in enumerate(pred_glob[:variant_count]):
            try:
                binder_seq = _sequence_from_structure_last_chain(cif_path)
                if binder_seq and len(binder_seq) >= 10:
                    p = _write_variant_or_fallback(i, binder_seq)
                    fastas.append(p)
            except Exception as ex:
                log.warning(f"Could not extract sequence from {cif_path}: {ex}")
        return fastas

    fastas = []
    for summary_path in summary_csvs:
        try:
            import pandas as pd
            df = pd.read_csv(summary_path)
            if "sequence" not in df.columns:
                seq_col = next((c for c in df.columns if "seq" in c.lower()), None)
                if not seq_col:
                    continue
                df = df.rename(columns={seq_col: "sequence"})
            for i, row in df.head(variant_count).iterrows():
                binder_seq = str(row.get("sequence", ""))
                if not binder_seq or len(binder_seq) < 10:
                    continue
                p = _write_variant_or_fallback(i, binder_seq)
                fastas.append(p)
        except Exception as ex:
            log.warning(f"Could not parse {summary_path}: {ex}")
    return fastas[:variant_count]
