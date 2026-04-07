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
import hashlib

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
) -> str:
    """Concise, generation-stable variant naming: Lxxxxxx_g01_v00"""
    tag = _short_lineage_tag(lineage_seed)
    return f"{tag}_g{generation_num:02d}_v{variant_index:02d}"


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


def _is_protein_chain(chain) -> bool:
    """Detect protein chain by backbone atoms (N, CA, C) or known residue names.
    Backbone-only CIFs may label all residues as UNK but still have protein backbone atoms."""
    _PROTEIN_BACKBONE = {"N", "CA", "C"}
    _RNA_ATOMS = {"O2'", "C2'", "C3'", "C4'", "O4'", "P", "OP1", "OP2"}
    for res in chain:
        if res.id[0] != " ":
            continue
        resname = res.get_resname().strip().upper()
        if _AA3_TO_1.get(resname):
            return True
        atom_names = {a.get_name().strip() for a in res}
        if _PROTEIN_BACKBONE.issubset(atom_names) and not _RNA_ATOMS.intersection(atom_names):
            return True
    return False


def _sequence_from_structure_last_chain(structure_path: str) -> str:
    """Extract sequence from the last protein chain (PXDesign outputs binder as final chain).
    Detects protein chains by backbone atoms, not just residue names — handles
    backbone-only CIFs where all residues are UNK."""
    from Bio.PDB import MMCIFParser, PDBParser
    ext = os.path.splitext(structure_path)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext == ".cif" else PDBParser(QUIET=True)
    struct = parser.get_structure("s", structure_path)

    protein_chain_id = None
    for chain in struct[0].get_chains():
        if _is_protein_chain(chain):
            protein_chain_id = chain.id

    if protein_chain_id is None:
        return ""
    return _sequence_from_structure(structure_path, protein_chain_id)


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
    """Check whether pxdesign CLI is reachable (returns True) or missing (returns False)."""
    import shutil
    pxdesign_bin = os.environ.get("PXDESIGN_CMD", "pxdesign")
    exe = pxdesign_bin.split()[0] if " " in pxdesign_bin else pxdesign_bin
    if os.path.isfile(exe):
        return True
    return shutil.which(exe) is not None


def _proteinmpnn_available() -> bool:
    """Check if standalone ProteinMPNN repo is available on this system."""
    mpnn_dir = os.environ.get("PROTEINMPNN_DIR", "")
    if mpnn_dir and os.path.isdir(mpnn_dir):
        return os.path.isfile(os.path.join(mpnn_dir, "protein_mpnn_run.py"))
    for candidate in ["/workspace/ProteinMPNN", "/opt/ProteinMPNN",
                      os.path.expanduser("~/ProteinMPNN")]:
        if os.path.isfile(os.path.join(candidate, "protein_mpnn_run.py")):
            os.environ["PROTEINMPNN_DIR"] = candidate
            return True
    return False


def _cif_to_enzyme_pdb(cif_path: str, output_pdb: str) -> bool:
    """
    Convert a PXDesign CIF to a PDB containing only the enzyme (protein) chain.

    PXDesign CIFs use auth_asym_id like 'A0', 'B0' etc. Biopython's MMCIFParser
    reads these as chain IDs. ProteinMPNN expects standard single-char PDB chain
    IDs ('A', 'B'). We also strip non-protein chains (RNA, ligands) since MPNN
    only designs protein sequences.
    """
    from Bio.PDB import MMCIFParser, PDBIO, Select

    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("s", cif_path)
    model = struct[0]

    chains = list(model.get_chains())
    if not chains:
        return False

    protein_chains = [c for c in chains if _is_protein_chain(c)]

    if not protein_chains:
        log.warning(f"No protein chains found in {cif_path}")
        return False

    remap = {}
    next_id = ord("A")
    for chain in protein_chains:
        old_id = chain.id
        new_id = chr(next_id)
        remap[old_id] = new_id
        next_id += 1
        if next_id > ord("Z"):
            break

    for chain in protein_chains:
        new_id = remap.get(chain.id)
        if new_id and new_id != chain.id:
            chain.id = new_id

    # ProteinMPNN requires standard amino acid names. Backbone-only CIFs from
    # PXDesign infer label every residue as UNK/XQB — replace with ALA so MPNN
    # can design real sequences for the backbone geometry.
    # CRITICAL: Biopython parses non-standard residues (UNK, XQB) as HETATM
    # with id[0] = "H_UNK". These get filtered out by PDBIO's default
    # accept_residue. We must fix the hetflag to " " (standard ATOM record).
    # Biopython uses res.id as a dict key in Chain, so we detach/re-add.
    n_renamed = 0
    for chain in protein_chains:
        residues_to_fix = []
        for res in chain:
            resname = res.get_resname().strip().upper()
            if resname not in _AA3_TO_1 or res.id[0] != " ":
                residues_to_fix.append(res)

        for res in residues_to_fix:
            old_id = res.id
            res.resname = "ALA"
            if old_id[0] != " ":
                new_id = (" ", old_id[1], old_id[2])
                if new_id in chain:
                    n_renamed += 1
                    continue
                chain.detach_child(old_id)
                res.id = new_id
                chain.add(res)
            n_renamed += 1
    if n_renamed:
        log.info(f"  Renamed {n_renamed} non-standard residues to ALA (ATOM) for ProteinMPNN")

    class ProteinOnlySelect(Select):
        def accept_chain(self, chain):
            return chain.id in remap.values()
        def accept_residue(self, residue):
            return residue.id[0] == " "

    io = PDBIO()
    io.set_structure(struct)
    io.save(output_pdb, ProteinOnlySelect())

    chain_info = ", ".join(f"{old}->{new}" for old, new in remap.items())
    # Verify PDB has ATOM records (not just HETATM)
    n_atom = 0
    if os.path.isfile(output_pdb):
        with open(output_pdb) as _pf:
            for _line in _pf:
                if _line.startswith("ATOM"):
                    n_atom += 1
    log.info(f"  CIF→PDB: extracted {len(protein_chains)} protein chain(s) [{chain_info}], {n_atom} ATOM records")
    if n_atom == 0:
        log.warning(f"  CIF→PDB produced 0 ATOM records — ProteinMPNN will fail!")
    return n_atom > 0


def _run_proteinmpnn_on_cif(cif_path: str, num_seqs: int = 4, temperature: float = 0.1) -> str:
    """
    Design a sequence for a backbone-only CIF using standalone ProteinMPNN.

    1. Convert CIF → PDB (enzyme-only, normalized chain IDs)
    2. Run ProteinMPNN's chain parser
    3. Run ProteinMPNN sequence design
    4. Return the best designed sequence (binder/last chain)
    """
    mpnn_dir = os.environ.get("PROTEINMPNN_DIR", "/workspace/ProteinMPNN")
    run_script = os.path.join(mpnn_dir, "protein_mpnn_run.py")

    import tempfile
    with tempfile.TemporaryDirectory(prefix="mpnn_") as tmp:
        # Step 1: CIF → enzyme-only PDB with clean chain IDs
        pdb_path = os.path.join(tmp, "enzyme.pdb")
        if not _cif_to_enzyme_pdb(cif_path, pdb_path):
            log.warning(f"CIF→PDB conversion failed for {cif_path}")
            return ""

        # Step 2: ProteinMPNN chain parsing
        input_dir = os.path.join(tmp, "pdbs")
        os.makedirs(input_dir, exist_ok=True)
        import shutil
        shutil.copy2(pdb_path, os.path.join(input_dir, "enzyme.pdb"))

        jsonl_path = os.path.join(tmp, "parsed_chains.jsonl")
        parse_script = os.path.join(mpnn_dir, "helper_scripts", "parse_multiple_chains.py")
        parse_cmd = [
            "python", parse_script,
            "--input_path", input_dir,
            "--output_path", jsonl_path,
        ]
        try:
            result = subprocess.run(parse_cmd, capture_output=True, text=True, timeout=120, check=True)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            log.warning(f"ProteinMPNN chain parsing failed: {e}")
            return ""
        except FileNotFoundError:
            log.warning(f"ProteinMPNN parse script not found: {parse_script}")
            return ""

        if not os.path.isfile(jsonl_path) or os.path.getsize(jsonl_path) == 0:
            log.warning("ProteinMPNN chain parser produced empty output")
            return ""

        with open(jsonl_path) as _jf:
            _parsed_content = _jf.read().strip()
        log.info(f"  MPNN parsed chains: {len(_parsed_content)} bytes, "
                 f"{_parsed_content.count(chr(10))+1} entries")

        # Step 3: run ProteinMPNN (enzyme-only, no RNA)
        out_dir = os.path.join(tmp, "output")
        os.makedirs(out_dir, exist_ok=True)
        mpnn_cmd = [
            "python", run_script,
            "--jsonl_path", jsonl_path,
            "--out_folder", out_dir,
            "--num_seq_per_target", str(num_seqs),
            "--sampling_temp", str(temperature),
            "--batch_size", str(min(8, num_seqs)),
        ]
        log.info(f"  MPNN cmd: {' '.join(mpnn_cmd)}")
        try:
            result = subprocess.run(mpnn_cmd, capture_output=True, text=True, timeout=300, check=True)
            if result.stdout:
                log.info(f"  MPNN stdout (tail): {result.stdout.strip()[-500:]}")
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            log.warning(f"ProteinMPNN sequence design failed: {e}")
            if hasattr(e, "stderr") and e.stderr:
                log.warning(f"  stderr: {e.stderr[-2000:]}")
            if hasattr(e, "stdout") and e.stdout:
                log.warning(f"  stdout: {e.stdout[-2000:]}")
            return ""
        except FileNotFoundError:
            log.warning(f"ProteinMPNN run script not found: {run_script}")
            return ""

        # Step 4: parse output — grab the last chain's sequence (binder)
        fa_glob = glob.glob(os.path.join(out_dir, "seqs", "*.fa"))
        if not fa_glob:
            fa_glob = glob.glob(os.path.join(out_dir, "**", "*.fa"), recursive=True)
        if not fa_glob:
            log.warning("ProteinMPNN produced no output FASTA files under %s", out_dir)
            return ""

        best_seq = ""
        best_score = float("inf")
        for fa_path in fa_glob:
            with open(fa_path) as f:
                lines = f.readlines()

            log.info(f"  MPNN FASTA {os.path.basename(fa_path)}: {len(lines)} lines, "
                     f"{sum(1 for l in lines if l.startswith('>'))} records")

            record_idx = 0
            i = 0
            while i < len(lines):
                header = lines[i].strip()
                i += 1
                seq_lines = []
                while i < len(lines) and not lines[i].startswith(">"):
                    seq_lines.append(lines[i].strip())
                    i += 1
                raw_seq = "".join(seq_lines)
                record_idx += 1
                if not raw_seq:
                    continue

                # Multi-chain outputs use '/' as chain separator — take
                # the last chain (binder) since we only want the designed
                # protein, not target/scaffold chains.
                if "/" in raw_seq:
                    chains_seqs = raw_seq.split("/")
                    seq = chains_seqs[-1]
                else:
                    seq = raw_seq

                x_frac = seq.count("X") / max(len(seq), 1)

                # Record 1 is always the input/native sequence echo.
                if record_idx == 1:
                    log.info(f"  MPNN rec 1 (input echo): {len(seq)}-aa, {x_frac:.0%} X | {header[:80]}")
                    continue

                score = float("inf")
                for part in header.split(","):
                    part = part.strip()
                    if part.startswith("score="):
                        try:
                            score = float(part.split("=")[1])
                        except ValueError:
                            pass

                log.info(f"  MPNN rec {record_idx}: {len(seq)}-aa, {x_frac:.0%} X, "
                         f"score={score:.3f} | {seq[:40]}...")

                if x_frac > 0.5:
                    log.info(f"  MPNN rec {record_idx}: >50%% X — skipping")
                    continue
                if score < best_score:
                    best_score = score
                    best_seq = seq

        if best_seq:
            log.info(f"  ProteinMPNN designed {len(best_seq)}-aa sequence (score={best_score:.2f})")
        else:
            log.warning("ProteinMPNN output FASTA contained no usable sequences")
        return best_seq


_MPNN_AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
_MPNN_AA_TO_IDX = {aa: i for i, aa in enumerate(_MPNN_AA_ORDER)}


def _build_pssm_jsonl(bias: dict, chain_a_len: int, chain_b_len: int,
                      rec_end: int, tmp_dir: str) -> str:
    """
    Convert RL bias dict ({"1-based_pos": {"AA": weight}}) into ProteinMPNN's
    PSSM JSONL format.

    ProteinMPNN expects: {pdb_name: {chain_id: {pssm_coef, pssm_bias, pssm_log_odds}}}
    where pssm_bias and pssm_log_odds are (L x 20) matrices.
    Chain A (scaffold) gets all zeros (no preference). Chain B (binder) gets RL bias.
    Amino acids in alphabetical order: ACDEFGHIKLMNPQRSTVWY.
    """
    import numpy as np
    pssm_a = np.zeros((chain_a_len, 20), dtype=float)
    pssm_b = np.zeros((chain_b_len, 20), dtype=float)

    for pos_str, aa_weights in bias.items():
        try:
            pos_0 = int(pos_str) - 1
        except ValueError:
            continue
        binder_idx = pos_0 - rec_end
        if 0 <= binder_idx < chain_b_len:
            for aa, weight in aa_weights.items():
                aa_idx = _MPNN_AA_TO_IDX.get(aa.upper())
                if aa_idx is not None:
                    pssm_b[binder_idx, aa_idx] = float(weight)

    coef_a = np.zeros(chain_a_len, dtype=float).tolist()
    coef_b = np.ones(chain_b_len, dtype=float).tolist()

    pssm_entry = {
        "backbone": {
            "A": {"pssm_coef": coef_a, "pssm_bias": pssm_a.tolist(),
                   "pssm_log_odds": pssm_a.tolist()},
            "B": {"pssm_coef": coef_b, "pssm_bias": pssm_b.tolist(),
                   "pssm_log_odds": pssm_b.tolist()},
        }
    }
    pssm_path = os.path.join(tmp_dir, "pssm.jsonl")
    with open(pssm_path, "w") as f:
        f.write(json.dumps(pssm_entry) + "\n")
    return pssm_path


def _run_mpnn_refinement(
    backbone_cif: str,
    current_seq: str,
    coords: dict,
    bias_json_path: str | None,
    variant_count: int,
    generation_num: int,
    temperature: float = 0.1,
) -> list[str]:
    """
    Iterative MPNN refinement: redesign a subset of linker positions on a fixed
    backbone, keeping HEPN/REC domains and most linker positions frozen.

    Each generation unfreezes more positions (exploration schedule) and applies
    RL bias weights as PSSM to guide which amino acids are preferred.

    Returns list of designed binder sequences (linker-only, for stitching).
    """
    if not _proteinmpnn_available():
        return []

    mpnn_dir = os.environ.get("PROTEINMPNN_DIR", "/workspace/ProteinMPNN")
    run_script = os.path.join(mpnn_dir, "protein_mpnn_run.py")

    import tempfile
    with tempfile.TemporaryDirectory(prefix="mpnn_refine_") as tmp:
        # Convert backbone CIF → PDB with ALA placeholders
        pdb_path = os.path.join(tmp, "backbone.pdb")
        if not _cif_to_enzyme_pdb(backbone_cif, pdb_path):
            log.warning("CIF→PDB conversion failed for refinement")
            return []

        # Parse chains
        input_dir = os.path.join(tmp, "pdbs")
        os.makedirs(input_dir)
        import shutil
        shutil.copy2(pdb_path, os.path.join(input_dir, "backbone.pdb"))

        jsonl_path = os.path.join(tmp, "parsed.jsonl")
        parse_script = os.path.join(mpnn_dir, "helper_scripts", "parse_multiple_chains.py")
        try:
            subprocess.run(
                ["python", parse_script, "--input_path", input_dir, "--output_path", jsonl_path],
                capture_output=True, text=True, timeout=120, check=True,
            )
        except Exception as e:
            log.warning(f"MPNN chain parsing failed for refinement: {e}")
            return []

        # Build fixed_positions JSONL: freeze everything except a subset of
        # linker positions. Unfreeze rate increases with generation.
        rec_end = coords["rec_end"]
        hepn1_start = coords["hepn1_start"]
        hepn1_end = coords["hepn1_end"]
        hepn2_start = coords["hepn2_start"]

        linker1_pos = list(range(rec_end, hepn1_start))
        linker2_pos = list(range(hepn1_end, hepn2_start))
        all_linker = linker1_pos + linker2_pos

        # Exploration schedule: unfreeze 10-50% of linker positions per gen
        import numpy as np
        rng = np.random.default_rng(generation_num * 7919)
        unfreeze_frac = min(0.10 + 0.05 * generation_num, 0.50)
        n_unfreeze = max(1, int(len(all_linker) * unfreeze_frac))

        # If RL bias is available, prefer unfreezing positions with strong signal
        unfreeze_weights = np.ones(len(all_linker), dtype=float)
        bias = {}
        if bias_json_path and os.path.exists(bias_json_path):
            try:
                with open(bias_json_path) as f:
                    bias = json.load(f)
                for idx, pos in enumerate(all_linker):
                    pos_str = str(pos + 1)
                    if pos_str in bias:
                        max_weight = max(bias[pos_str].values())
                        unfreeze_weights[idx] = 1.0 + abs(max_weight)
            except Exception:
                pass

        unfreeze_weights = np.clip(unfreeze_weights, 0.01, None)
        unfreeze_weights /= unfreeze_weights.sum()
        unfreeze_indices = rng.choice(len(all_linker), size=min(n_unfreeze, len(all_linker)),
                                      replace=False, p=unfreeze_weights)
        designable_positions = {all_linker[i] for i in unfreeze_indices}

        # MPNN fixed_positions_jsonl: {pdb_name: {chain: [1-based int positions]}}
        # Fixed = all positions NOT in designable set. Positions are 1-based
        # relative to the chain as parsed by MPNN.
        from Bio.PDB import PDBParser as _PDBParser
        _pdb_struct = _PDBParser(QUIET=True).get_structure("b", pdb_path)
        _chains = {c.id: sum(1 for r in c.get_residues() if r.id[0] == ' ')
                   for c in _pdb_struct[0].get_chains()}
        binder_len = _chains.get("B", coords["binder_length"])
        chain_a_len = _chains.get("A", 0)

        fixed_in_A = list(range(1, chain_a_len + 1))
        fixed_in_B = [p - rec_end + 1 for p in range(rec_end, rec_end + binder_len)
                      if p not in designable_positions]
        fixed_positions = {"backbone": {"A": fixed_in_A, "B": fixed_in_B}}
        fixed_pos_path = os.path.join(tmp, "fixed_positions.jsonl")
        with open(fixed_pos_path, "w") as f:
            f.write(json.dumps(fixed_positions) + "\n")

        log.info(f"  [MPNN refine] gen={generation_num}: unfreezing {len(designable_positions)}/{len(all_linker)} "
                 f"linker positions ({unfreeze_frac:.0%}), {len(bias)} RL bias entries")

        out_dir = os.path.join(tmp, "output")
        os.makedirs(out_dir)

        design_temp = max(0.1, 0.3 - 0.02 * generation_num)

        # Build PSSM bias from RL weights for MPNN-native soft guidance
        pssm_path = None
        if bias:
            pssm_path = _build_pssm_jsonl(bias, chain_a_len, binder_len, rec_end, tmp)

        mpnn_cmd = [
            "python", run_script,
            "--jsonl_path", jsonl_path,
            "--out_folder", out_dir,
            "--num_seq_per_target", str(variant_count),
            "--sampling_temp", str(design_temp),
            "--batch_size", str(min(8, variant_count)),
            "--fixed_positions_jsonl", fixed_pos_path,
        ]
        if pssm_path:
            mpnn_cmd.extend(["--pssm_jsonl", pssm_path, "--pssm_multi", "0.5"])

        try:
            result = subprocess.run(mpnn_cmd, capture_output=True, text=True, timeout=300, check=True)
        except subprocess.CalledProcessError as e:
            stderr_tail = (e.stderr or "")[-1500:]
            log.warning(f"MPNN refinement failed (exit {e.returncode}): {stderr_tail}")
            return []
        except Exception as e:
            log.warning(f"MPNN refinement failed: {e}")
            return []

        # Parse output FASTAs — same logic as _run_proteinmpnn_on_cif
        fa_glob = glob.glob(os.path.join(out_dir, "seqs", "*.fa"))
        if not fa_glob:
            return []

        designed_seqs = []
        for fa_path in fa_glob:
            with open(fa_path) as f:
                lines = f.readlines()
            record_idx = 0
            i = 0
            while i < len(lines):
                header = lines[i].strip()
                i += 1
                seq_lines = []
                while i < len(lines) and not lines[i].startswith(">"):
                    seq_lines.append(lines[i].strip())
                    i += 1
                raw_seq = "".join(seq_lines)
                record_idx += 1
                if not raw_seq or record_idx == 1:
                    continue
                seq = raw_seq.split("/")[-1] if "/" in raw_seq else raw_seq
                if seq.count("X") / max(len(seq), 1) > 0.5:
                    continue
                designed_seqs.append(seq)

        log.info(f"  [MPNN refine] Generated {len(designed_seqs)} refined sequences")
        return designed_seqs[:variant_count]


def _apply_bias_to_sequence(full_seq: str, bias_json_path: str, coords: dict) -> str:
    """
    Apply RL bias matrix to the designed full sequence.
    Only modifies linker regions (preserves HEPN catalytic domains and REC).
    Bias matrix format: {"position_1based": {"amino_acid": weight, ...}, ...}
    Higher positive weight → prefer this AA; substitution applied when weight > threshold.
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
            pos = int(pos_str) - 1  # Mutations are 1-based; convert to 0-based index
        except ValueError:
            continue
        if pos < 0 or pos >= len(seq_list):
            continue
        # Only modify linker regions (between REC↔HEPN1 and HEPN1↔HEPN2)
        in_linker1 = rec_end <= pos < hepn1_start
        in_linker2 = hepn1_end <= pos < hepn2_start
        if not (in_linker1 or in_linker2):
            continue
        # Pick the amino acid with highest positive bias weight
        best_aa = max(aa_weights, key=aa_weights.get)
        best_weight = aa_weights[best_aa]
        if best_weight > 0.5 and best_aa != seq_list[pos]:
            seq_list[pos] = best_aa
            mutations_applied += 1

    if mutations_applied > 0:
        log.info(f"  RL bias applied {mutations_applied} substitution(s) in linker regions")
    return "".join(seq_list)


def _sequence_from_fasta_or_json(fasta_path: str, json_path: str, variant_id: str) -> str:
    """Get baseline sequence from FASTA or base JSON.

    Searches multiple candidate JSON paths in case the base_json_dir doesn't
    align with the actual jsons/ directory on this machine.
    """
    if fasta_path and os.path.exists(fasta_path):
        with open(fasta_path) as f:
            return "".join(l.strip() for l in f if not l.startswith(">"))

    candidates = [json_path]
    scripts_dir = os.path.dirname(os.path.abspath(__file__))
    candidates.append(os.path.join(scripts_dir, "..", "jsons", f"{variant_id}.json"))
    candidates.append(os.path.join(scripts_dir, "..", "outputs", "jsons", f"{variant_id}.json"))

    for path in candidates:
        if not path or not os.path.exists(path):
            continue
        try:
            with open(path) as f:
                data = json.load(f)
            ent = data[0]["sequences"][0]
            prot = ent.get("proteinChain", ent.get("protein", {}))
            seq = prot.get("sequence", "")
            if seq:
                return seq
        except (json.JSONDecodeError, KeyError, IndexError):
            continue
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
    # All coords are 0-based Python-slice indices; linker spans [rec_end, hepn1_start)
    linker1_len = max(0, hepn1_start - rec_end)
    linker2_len = max(0, hepn2_start - hepn1_end)
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
    Runs PXDesign via: pxdesign infer -i <yaml> -o <dir> --N_sample N
    Generation only (no evaluation); we compute scores ourselves with Protenix later.
    Designs linkers only; stitches wild-type HEPN1/HEPN2 into output.
    Raises RuntimeError if PXDesign is unavailable or fails.
    """
    log.info(f"Generating {variant_count} designs for {variant_id}...")
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
    if not full_seq:
        log.warning("Could not read baseline sequence from any source")
    seq_deficit = coords["hepn2_end"] - len(full_seq)
    if seq_deficit > 0:
        log.warning(f"Baseline sequence ({len(full_seq)} aa) shorter than HEPN2 end ({coords['hepn2_end']}); "
                     f"padding {seq_deficit} residues with Glycine")
        full_seq += "G" * seq_deficit

    abs_out = os.path.abspath(output_dir)
    name_seed = lineage_seed or variant_id

    def _write_variant_fasta(i: int, variant_full_seq: str) -> str:
        name = _compact_variant_name(
            lineage_seed=name_seed,
            generation_num=max(1, int(generation_num)),
            variant_index=i,
        )
        fasta_path = os.path.join(output_dir, f"{name}.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{name}\n{variant_full_seq}\n")
        return fasta_path

    # ── Generation strategy: PXDesign for gen 1 backbone, MPNN refinement for gen 2+ ──
    backbone_cache_dir = os.path.join(os.path.dirname(abs_out), "backbone_cache")
    os.makedirs(backbone_cache_dir, exist_ok=True)
    backbone_cif_cache = os.path.join(backbone_cache_dir, f"{name_seed}_backbone.cif")

    if generation_num <= 1 or not os.path.isfile(backbone_cif_cache):
        # Gen 1: run PXDesign to generate novel backbone geometry
        if not _pxdesign_available():
            raise RuntimeError(
                "PXDesign CLI not found. Install PXDesign (bash install.sh --env pxdesign) "
                "and set PXDESIGN_CMD to the binary path. See setup_dual_env.sh."
            )

        succeeded = _run_pxdesign_cli(
            baseline_structure, variant_id, coords, output_dir, variant_count, bias_json_path,
        )
        if not succeeded:
            raise RuntimeError(
                "PXDesign exited with non-zero status. Check logs above for stderr. "
                "Common fix: pip install deepspeed==0.14.5 in the pxdesign env."
            )

        fastas = _parse_pxdesign_outputs(
            abs_out, variant_count, full_seq, coords, bias_json_path,
            output_dir, name_seed, generation_num,
        )
        if not fastas:
            raise RuntimeError(
                f"PXDesign ran successfully but produced no usable variants under {abs_out}. "
                "Check design_outputs/*/summary.csv or predictions/*.cif."
            )

        # Cache the best backbone CIF for future refinement generations
        pred_cifs = glob.glob(os.path.join(abs_out, "**", "predictions", "*.cif"), recursive=True)
        if pred_cifs:
            import shutil
            shutil.copy2(pred_cifs[0], backbone_cif_cache)
            log.info(f"  Cached backbone CIF for future MPNN refinement: {backbone_cif_cache}")

        return fastas
    else:
        # Gen 2+: MPNN refinement on the cached PXDesign backbone
        log.info(f"  [Gen {generation_num}] MPNN iterative refinement on cached backbone")
        refined_seqs = _run_mpnn_refinement(
            backbone_cif=backbone_cif_cache,
            current_seq=full_seq,
            coords=coords,
            bias_json_path=bias_json_path,
            variant_count=variant_count,
            generation_num=generation_num,
        )

        if not refined_seqs:
            log.warning("MPNN refinement produced no sequences; falling back to PXDesign")
            succeeded = _run_pxdesign_cli(
                baseline_structure, variant_id, coords, output_dir, variant_count, bias_json_path,
            )
            if succeeded:
                return _parse_pxdesign_outputs(
                    abs_out, variant_count, full_seq, coords, bias_json_path,
                    output_dir, name_seed, generation_num,
                )
            raise RuntimeError("Both MPNN refinement and PXDesign failed.")

        from utils.hepn_structural_stitch import stitch_hepn_into_binder
        fastas = []
        for i, binder_seq in enumerate(refined_seqs):
            binder_seq = _resolve_unknown_residues(binder_seq, full_seq[coords["rec_end"]:])
            full = stitch_hepn_into_binder(binder_seq, full_seq, coords)
            if not full:
                continue
            full = _resolve_unknown_residues(full, full_seq)
            full = _apply_bias_to_sequence(full, bias_json_path, coords)
            name = _compact_variant_name(name_seed, max(1, int(generation_num)), i)
            fasta_path = os.path.join(output_dir, f"{name}.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">{name}\n{full}\n")
            fastas.append(fasta_path)

        if not fastas:
            raise RuntimeError("MPNN refinement produced sequences but stitching failed for all.")

        log.info(f"  [MPNN refine] Wrote {len(fastas)} variant FASTA(s)")
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
            log.error(
                "PXDesign FAILED (exit %d). Troubleshooting:\n"
                "  1. Try PXDESIGN_SUBCOMMAND=infer if your version doesn't support 'pipeline'\n"
                "  2. Check GPU memory: PXDesign + Protenix may exceed VRAM\n"
                "  3. Verify: PXDESIGN_CMD=%s  PXDESIGN_SUBCOMMAND=%s\n"
                "  4. Run manually: %s",
                result.returncode, pxdesign_bin, sub, " ".join(cmd),
            )
            if err:
                log.error("PXDesign stderr (tail):\n%s", err[-tail:] if len(err) > tail else err)
            if out:
                log.error("PXDesign stdout (tail):\n%s", out[-tail:] if len(out) > tail else out)
            return False
        log.info("PXDesign completed successfully (exit 0)")
        return True
    except subprocess.TimeoutExpired:
        log.error("PXDesign timed out after 1 hour — consider reducing --N_sample or using 'preview' preset")
        return False
    except FileNotFoundError:
        log.error(
            "PXDesign binary not found: %s\n"
            "  Install PXDesign and set PXDESIGN_CMD to the binary path.\n"
            "  Or run: conda activate pxdesign && which pxdesign",
            pxdesign_exec[0],
        )
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

    def _write_stitched_variant(i: int, binder_seq: str) -> str | None:
        binder_seq = _resolve_unknown_residues(binder_seq, full_seq[coords["rec_end"]:])
        full = stitch_hepn_into_binder(binder_seq, full_seq, coords)
        if not full:
            log.warning(f"Stitching failed for design {i}; skipping (no fallback)")
            return None
        full = _resolve_unknown_residues(full, full_seq)
        full = _apply_bias_to_sequence(full, bias_json_path, coords)
        name = _compact_variant_name(name_seed, max(1, int(generation_num)), i)
        fasta_path = os.path.join(output_dir, f"{name}.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{name}\n{full}\n")
        return fasta_path

    sample_csvs = glob.glob(os.path.join(abs_out, "**", "sample_level_output.csv"), recursive=True)
    design_out = os.path.join(abs_out, "design_outputs")
    summary_csvs = sample_csvs or glob.glob(os.path.join(design_out, "**", "summary.csv"), recursive=True)

    if not summary_csvs:
        log.warning("No sequence CSV found; checking CIF predictions")
        pred_glob = glob.glob(os.path.join(abs_out, "**", "predictions", "*.cif"), recursive=True)
        if not pred_glob:
            return []

        mpnn_available = _proteinmpnn_available()
        fastas = []
        for i, cif_path in enumerate(pred_glob[:variant_count]):
            try:
                binder_seq = _sequence_from_structure_last_chain(cif_path)
                if not binder_seq or len(binder_seq) < 10:
                    continue
                x_frac = binder_seq.count("X") / len(binder_seq)
                if x_frac > 0.5:
                    if mpnn_available:
                        log.info(f"CIF {cif_path}: backbone-only ({x_frac:.0%} X) — running ProteinMPNN for sequence design")
                        mpnn_seq = _run_proteinmpnn_on_cif(cif_path)
                        if mpnn_seq:
                            mpnn_x = mpnn_seq.count("X")
                            log.info(f"  MPNN returned {len(mpnn_seq)}-aa, {mpnn_x} X residues ({mpnn_seq[:30]}...)")
                            binder_seq = mpnn_seq
                        else:
                            log.warning(f"ProteinMPNN failed on {cif_path}; resolving from baseline")
                            binder_seq = _sequence_from_structure_last_chain(cif_path)
                    else:
                        log.info(
                            f"CIF {cif_path}: backbone-only ({x_frac:.0%} X). "
                            "Resolving from baseline (install ProteinMPNN for better sequence design)."
                        )
                path = _write_stitched_variant(i, binder_seq)
                if path:
                    fastas.append(path)
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
                path = _write_stitched_variant(i, binder_seq)
                if path:
                    fastas.append(path)
        except Exception as ex:
            log.warning(f"Could not parse {summary_path}: {ex}")
    return fastas[:variant_count]
