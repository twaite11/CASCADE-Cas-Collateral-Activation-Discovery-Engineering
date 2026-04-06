"""
PXDesign wrapper for CASCADE evolution.
Uses the actual PXDesign CLI: pxdesign infer -i <yaml> -o <dir> --N_sample N
Generation only; we compute Protenix scores downstream. Generates YAML from baseline structure + metadata,
runs inference, parses CIF output to variant FASTAs.
Falls back to built-in sequence-level linker mutator when PXDesign CLI is unavailable.
"""
import subprocess
import json
import os
import glob
import logging
import hashlib
from pathlib import Path

log = logging.getLogger(__name__)

# Map 3-letter to 1-letter amino acid codes (standard + common)
_AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _short_lineage_tag(seed: str, fallback: str = "lineage") -> str:
    """Compact deterministic lineage tag."""
    seed = (seed or "").strip()
    if not seed:
        seed = fallback
    digest = hashlib.sha1(seed.encode("utf-8")).hexdigest()[:6]
    return f"L{digest}"


def _compact_variant_name(
    lineage_seed: str,
    generation_num: int,
    variant_index: int,
    fallback: bool = False,
) -> str:
    tag = _short_lineage_tag(lineage_seed)
    base = f"{tag}_g{generation_num:02d}_v{variant_index:02d}"
    if fallback:
        return f"{base}_fb"
    return base


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
    ext = os.path.splitext(structure_path)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext == ".cif" else PDBParser(QUIET=True)
    struct = parser.get_structure("s", structure_path)
    chain = struct[0][chain_id]
    seq = []
    for res in chain:
        if res.id[0] != " ":
            continue
        resname = res.get_resname().strip().upper()
        seq.append(_AA3_TO_1.get(resname, "X"))
    return "".join(seq)


def _sequence_from_structure_last_chain(structure_path: str) -> str:
    """Extract sequence from the last chain (PXDesign outputs binder as final chain)."""
    chain_ids = _get_structure_chain_ids(structure_path)
    if not chain_ids:
        return ""
    return _sequence_from_structure(structure_path, chain_ids[-1])


def _resolve_unknown_residues(seq: str, baseline_seq: str = "") -> str:
    """
    Replace X (unknown) residues in designed sequences with valid amino acids.
    Strategy: use the baseline residue at that position if available, otherwise
    substitute Glycine (smallest, least disruptive to backbone geometry).
    Returns the resolved sequence with no X characters.
    """
    if "X" not in seq:
        return seq
    resolved = []
    for i, aa in enumerate(seq):
        if aa == "X":
            if i < len(baseline_seq) and baseline_seq[i] not in ("X", "-", ""):
                resolved.append(baseline_seq[i])
            else:
                resolved.append("G")  # Glycine fallback
        else:
            resolved.append(aa)
    result = "".join(resolved)
    n_resolved = seq.count("X")
    if n_resolved > 0:
        log.info(f"  Resolved {n_resolved} unknown (X) residue(s) in designed sequence")
    return result


def _pxdesign_available() -> bool:
    """Check whether pxdesign CLI is reachable."""
    import shutil
    pxdesign_bin = os.environ.get("PXDESIGN_CMD", "pxdesign")
    exe = pxdesign_bin.split()[0] if " " in pxdesign_bin else pxdesign_bin
    if os.path.isfile(exe):
        return True
    return shutil.which(exe) is not None


_AA_ALL_STANDARD = list("ACDEFGHIKLMNPQRSTVWY")
_LINKER_POOL = list("AGSTDENQKR")


def _generate_variants_sequence_mutator(
    full_seq: str,
    coords: dict,
    variant_count: int,
    bias_json_path: str | None,
    generation_num: int,
    seed: int = 42,
) -> list[str]:
    """
    Pure-sequence linker mutator: produces variant full-length sequences by mutating
    only linker residues (REC->HEPN1 gap and HEPN1->HEPN2 gap), preserving catalytic
    domains exactly. Uses RL bias when available; otherwise applies conservative
    substitutions from a linker-friendly amino acid pool.

    Returns list of mutated full-length sequences (same length as full_seq).
    """
    import numpy as np
    rng = np.random.default_rng(seed + generation_num)

    rec_end = coords["rec_end"]
    hepn1_start = coords["hepn1_start"]
    hepn1_end = coords["hepn1_end"]
    hepn2_start = coords["hepn2_start"]

    linker1_positions = list(range(rec_end, hepn1_start))
    linker2_positions = list(range(hepn1_end, hepn2_start))
    mutable_positions = linker1_positions + linker2_positions

    if not mutable_positions:
        log.warning("No mutable linker positions — domains are contiguous. Returning baseline.")
        return [full_seq] * variant_count

    bias = {}
    if bias_json_path and os.path.exists(bias_json_path):
        try:
            with open(bias_json_path) as f:
                bias = json.load(f)
            log.info(f"  Sequence mutator: loaded RL bias ({len(bias)} positions)")
        except Exception:
            pass

    n_mutable = len(mutable_positions)
    base_rate = min(0.05 + 0.01 * generation_num, 0.20)

    variants = []
    for vi in range(variant_count):
        seq_list = list(full_seq)

        n_mutations = max(1, rng.poisson(base_rate * n_mutable))
        n_mutations = min(n_mutations, n_mutable)

        chosen = rng.choice(mutable_positions, size=n_mutations, replace=False)

        for pos in chosen:
            pos_str = str(pos + 1)
            current_aa = seq_list[pos]

            if pos_str in bias:
                aa_weights = bias[pos_str]
                candidates = [(aa, w) for aa, w in aa_weights.items()
                              if aa != current_aa and w > 0 and aa in _AA_ALL_STANDARD]
                if candidates:
                    aas, ws = zip(*candidates)
                    ws_arr = np.array(ws, dtype=float)
                    ws_arr /= ws_arr.sum()
                    seq_list[pos] = rng.choice(list(aas), p=ws_arr)
                    continue

            pool = [aa for aa in _LINKER_POOL if aa != current_aa]
            seq_list[pos] = rng.choice(pool) if pool else current_aa

        variants.append("".join(seq_list))

    return variants


def _apply_bias_to_sequence(full_seq: str, bias_json_path: str, coords: dict) -> str:
    """
    Apply RL bias matrix to the designed full sequence.
    Only modifies linker regions (preserves HEPN catalytic domains and REC).
    Bias matrix format: {"position_1based": {"amino_acid": weight, ...}, ...}
    Higher positive weight -> prefer this AA; substitution applied when weight > threshold.
    """
    if not bias_json_path or not os.path.exists(bias_json_path):
        return full_seq
    with open(bias_json_path) as f:
        bias = json.load(f)
    if not bias:
        return full_seq

    rec_end = coords["rec_end"]
    hepn1_start = coords["hepn1_start"]
    hepn1_end = coords["hepn1_end"]
    hepn2_start = coords["hepn2_start"]

    seq_list = list(full_seq)
    mutations_applied = 0

    for pos_str, aa_weights in bias.items():
        try:
            pos = int(pos_str) - 1
        except ValueError:
            continue
        if pos < 0 or pos >= len(seq_list):
            continue
        in_linker1 = rec_end <= pos < hepn1_start
        in_linker2 = hepn1_end <= pos < hepn2_start
        if not (in_linker1 or in_linker2):
            continue
        best_aa = max(aa_weights, key=aa_weights.get)
        best_weight = aa_weights[best_aa]
        if best_weight > 0.5 and best_aa != seq_list[pos]:
            seq_list[pos] = best_aa
            mutations_applied += 1

    if mutations_applied > 0:
        log.info(f"  RL bias applied {mutations_applied} substitution(s) in linker regions")
    return "".join(seq_list)


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

    linker1_len = max(0, hepn1_start - rec_end)
    linker2_len = max(0, hepn2_start - hepn1_end)
    binder_length = linker1_len + linker2_len
    MIN_BINDER_LENGTH = 20
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
    """Build PXDesign YAML: protein-only target (REC crop), linker-only binder_length."""
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
    generation_num: int = 1,
    lineage_seed: str | None = None,
):
    """
    Runs PXDesign via: pxdesign pipeline -i <yaml> -o <dir> --N_sample N
    Generation only (no evaluation); we compute scores ourselves with Protenix later.
    Designs linkers only; stitches wild-type HEPN1/HEPN2 into output.
    Falls back to built-in sequence-level linker mutator when PXDesign CLI is unavailable.
    """
    log.info(f"Generating {variant_count} designs for {variant_id}...")
    os.makedirs(output_dir, exist_ok=True)

    coords = generate_frozen_rec_config(metadata_path, variant_id, metadata_override)

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

    abs_out = os.path.abspath(output_dir)
    name_seed = lineage_seed or variant_id

    def _write_variant_fasta(i: int, variant_full_seq: str, is_fallback: bool = False) -> str:
        name = _compact_variant_name(
            lineage_seed=name_seed,
            generation_num=max(1, int(generation_num)),
            variant_index=i,
            fallback=is_fallback,
        )
        fasta_path = os.path.join(output_dir, f"{name}.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{name}\n{variant_full_seq}\n")
        return fasta_path

    # ── Decide engine: PXDesign CLI or built-in sequence mutator ──
    use_pxdesign = _pxdesign_available() and os.environ.get("CASCADE_FORCE_SEQ_MUTATOR", "") != "1"
    pxdesign_succeeded = False

    if use_pxdesign:
        pxdesign_succeeded = _run_pxdesign_cli(
            baseline_structure, variant_id, coords, output_dir, variant_count, bias_json_path,
        )

    if use_pxdesign and pxdesign_succeeded:
        fastas = _parse_pxdesign_outputs(
            abs_out, variant_count, full_seq, coords, bias_json_path,
            output_dir, name_seed, generation_num,
        )
        if fastas:
            return fastas
        log.warning("PXDesign ran but produced no usable variants; falling through to sequence mutator.")

    # ── Sequence-level linker mutator (no PXDesign dependency) ──
    if not use_pxdesign:
        log.info(
            "[SeqMutator] PXDesign CLI not found — using built-in sequence-level linker mutator. "
            "Set PXDESIGN_CMD to enable structure-guided generation."
        )
    else:
        log.info("[SeqMutator] Falling back to sequence-level linker mutator after PXDesign failure.")

    seed_val = hash(variant_id) & 0xFFFFFFFF
    mutant_seqs = _generate_variants_sequence_mutator(
        full_seq, coords, variant_count, bias_json_path, generation_num, seed=seed_val,
    )

    fastas = []
    for i, mseq in enumerate(mutant_seqs):
        mseq = _resolve_unknown_residues(mseq, full_seq)
        mseq = _apply_bias_to_sequence(mseq, bias_json_path, coords)
        fastas.append(_write_variant_fasta(i, mseq, is_fallback=False))
    log.info(f"[SeqMutator] Wrote {len(fastas)} variant FASTA(s) to {output_dir}")
    return fastas


def _run_pxdesign_cli(
    baseline_structure: str,
    variant_id: str,
    coords: dict,
    output_dir: str,
    variant_count: int,
    bias_json_path: str | None,
) -> bool:
    """Invoke the PXDesign CLI subprocess. Returns True if the process exited 0."""
    yaml_path = _build_pxdesign_yaml(baseline_structure, variant_id, coords, output_dir)
    abs_out = os.path.abspath(output_dir)

    pxdesign_bin = os.environ.get("PXDESIGN_CMD", "pxdesign")
    pxdesign_exec = pxdesign_bin.split() if " " in pxdesign_bin else [pxdesign_bin]

    sub = os.environ.get("PXDESIGN_SUBCOMMAND", "pipeline").strip().lower()
    if sub not in ("pipeline", "infer"):
        log.warning(f"PXDESIGN_SUBCOMMAND={sub!r} unknown; using 'pipeline'")
        sub = "pipeline"

    cmd = [
        *pxdesign_exec, sub,
        "-i", yaml_path,
        "-o", abs_out,
        "--N_sample", str(variant_count),
        "--dtype", "bf16",
    ]
    if sub == "pipeline":
        cmd.extend(["--preset", "preview", "--N_max_runs", "1"])
    if bias_json_path and os.path.exists(bias_json_path):
        log.info(f"  -> RL bias matrix loaded; will apply to designed sequences post-stitching")

    log.info(f"  -> PXDesign cmd: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,
            cwd=os.path.dirname(os.path.abspath(__file__)) or ".",
        )
        if result.returncode != 0:
            tail = 6000
            err = (result.stderr or "").strip()
            out = (result.stdout or "").strip()
            log.warning(f"PXDesign exited with code {result.returncode}")
            if err:
                log.error("PXDesign stderr (tail):\n%s", err[-tail:] if len(err) > tail else err)
            if out:
                log.error("PXDesign stdout (tail):\n%s", out[-tail:] if len(out) > tail else out)
            return False
        return True
    except subprocess.TimeoutExpired:
        log.error("PXDesign timed out after 1 hour")
        return False
    except FileNotFoundError:
        log.error("PXDesign binary not found: %s", pxdesign_exec[0])
        return False


def _parse_pxdesign_outputs(
    abs_out: str,
    variant_count: int,
    full_seq: str,
    coords: dict,
    bias_json_path: str | None,
    output_dir: str,
    name_seed: str,
    generation_num: int,
) -> list[str]:
    """Parse PXDesign CSV/CIF outputs, stitch HEPN domains, write variant FASTAs."""
    from utils.hepn_structural_stitch import stitch_hepn_into_binder

    def _write_variant_or_fallback(i: int, binder_seq: str) -> str:
        binder_seq = _resolve_unknown_residues(binder_seq, full_seq[coords["rec_end"]:])
        full = stitch_hepn_into_binder(binder_seq, full_seq, coords)
        if full:
            full = _resolve_unknown_residues(full, full_seq)
            full = _apply_bias_to_sequence(full, bias_json_path, coords)
            name = _compact_variant_name(name_seed, max(1, int(generation_num)), i, fallback=False)
            fasta_path = os.path.join(output_dir, f"{name}.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">{name}\n{full}\n")
            return fasta_path
        fallback_name = _compact_variant_name(name_seed, max(1, int(generation_num)), i, fallback=True)
        fallback_path = os.path.join(output_dir, f"{fallback_name}.fasta")
        with open(fallback_path, "w") as f:
            f.write(f">{fallback_name}\n{full_seq}\n")
        log.warning(f"Stitching failed for design {i}; writing baseline as fallback (will receive penalty)")
        return fallback_path

    sample_csvs = glob.glob(os.path.join(abs_out, "**", "sample_level_output.csv"), recursive=True)
    design_out = os.path.join(abs_out, "design_outputs")
    summary_csvs = sample_csvs or glob.glob(os.path.join(design_out, "**", "summary.csv"), recursive=True)

    if not summary_csvs:
        log.warning("No sequence CSV found; checking CIF predictions")
        pred_glob = glob.glob(os.path.join(abs_out, "**", "predictions", "*.cif"), recursive=True)
        if not pred_glob:
            return []
        fastas = []
        for i, cif_path in enumerate(pred_glob[:variant_count]):
            try:
                binder_seq = _sequence_from_structure_last_chain(cif_path)
                if not binder_seq or len(binder_seq) < 10:
                    continue
                x_frac = binder_seq.count("X") / len(binder_seq)
                if x_frac > 0.5:
                    log.error(
                        f"CIF {cif_path}: {x_frac:.0%} of binder residues are unknown (X). "
                        "Backbone-only outputs have no sequence identity."
                    )
                    continue
                fastas.append(_write_variant_or_fallback(i, binder_seq))
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
                fastas.append(_write_variant_or_fallback(i, binder_seq))
        except Exception as ex:
            log.warning(f"Could not parse {summary_path}: {ex}")
    return fastas[:variant_count]
