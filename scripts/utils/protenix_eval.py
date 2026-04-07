import os
import json
import subprocess
import glob
import random
import shutil
import logging

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Eval engine detection: cattle-prod (Rust) preferred, protenix (Python) fallback
# ---------------------------------------------------------------------------
# Set EVAL_CMD to override auto-detection (e.g. EVAL_CMD=protenix or EVAL_CMD=/path/to/cattle-prod)
_EVAL_CMD_OVERRIDE = os.environ.get("EVAL_CMD", "").strip()

def _detect_eval_engine():
    """Detect the best available structure prediction engine.
    Returns (binary, engine_name) where engine_name is 'cattle-prod' or 'protenix'."""
    if _EVAL_CMD_OVERRIDE:
        name = "cattle-prod" if "cattle-prod" in _EVAL_CMD_OVERRIDE else "protenix"
        return _EVAL_CMD_OVERRIDE, name

    cp_bin = shutil.which("cattle-prod")
    if not cp_bin:
        for candidate in [
            os.path.join(os.path.dirname(__file__), "..", "..", "..", "model_rustprot", "cattle-prod", "target", "release", "cattle-prod"),
            os.path.join(os.path.dirname(__file__), "..", "rust", "cattle-prod", "target", "release", "cattle-prod"),
        ]:
            if os.name == "nt":
                candidate += ".exe"
            if os.path.isfile(candidate):
                cp_bin = os.path.abspath(candidate)
                break

    if cp_bin:
        return cp_bin, "cattle-prod"

    return "protenix", "protenix"

EVAL_BIN, EVAL_ENGINE = _detect_eval_engine()
log.info(f"CASCADE eval engine: {EVAL_ENGINE} ({EVAL_BIN})")
_CATTLE_PROD_STRICT = os.environ.get("CATTLE_PROD_STRICT", "1").strip().lower() not in {"0", "false", "no"}

# RNA Constants
# Default dummy spacer for backward compatibility. Prefer loading real
# fusion targets from data/fusion_targets.json via load_fusion_target().
_RNA_COMPLEMENT = str.maketrans("AUGC", "UACG")

_FUSION_TARGETS_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "data", "fusion_targets.json")
_loaded_target = None


def _spacer_to_target(spacer_rna: str) -> str:
    """Convert a 24nt spacer RNA to its target RNA (reverse complement + flanks)."""
    rc = spacer_rna[::-1].translate(_RNA_COMPLEMENT)
    return "AAAAAA" + rc + "AAAAAA"


def load_fusion_target(target_id: str = None) -> tuple:
    """Load a fusion target from data/fusion_targets.json.
    Returns (spacer_rna, target_rna, target_info_dict).
    Falls back to dummy if file missing or target not found."""
    global _loaded_target
    if _loaded_target and (target_id is None or _loaded_target[2].get("id") == target_id):
        return _loaded_target

    if os.path.exists(_FUSION_TARGETS_PATH):
        try:
            with open(_FUSION_TARGETS_PATH) as f:
                data = json.load(f)
            targets = {t["id"]: t for t in data.get("targets", [])}
            tid = target_id or os.environ.get("CASCADE_FUSION_TARGET") or data.get("default_target_id")
            if tid and tid in targets:
                info = targets[tid]
                spacer = info["spacer_rna"]
                target = _spacer_to_target(spacer)
                _loaded_target = (spacer, target, info)
                log.info(f"Loaded fusion target: {info['fusion']} ({tid}) — {info['cancer']}")
                return _loaded_target
        except Exception as e:
            log.warning(f"Could not load fusion targets: {e}")

    spacer = "GUCGACUGACGUACGUACGUACGU"
    target = _spacer_to_target(spacer)
    _loaded_target = (spacer, target, {"id": "dummy", "fusion": "dummy", "cancer": "N/A"})
    return _loaded_target


DUMMY_SPACER_RNA, DUMMY_TARGET_RNA, _TARGET_INFO = load_fusion_target()
TARGET_REGION = DUMMY_SPACER_RNA[::-1].translate(_RNA_COMPLEMENT)

# Optimal spacer lengths per Cas13 subtype (from literature)
_SUBTYPE_SPACER_LEN = {
    "cas13a": 28,
    "cas13b": 30,
    "cas13d": 23,
}
_DEFAULT_SPACER_LEN = 24


def get_spacer_for_subtype(subtype: str = "unknown") -> str:
    """Return the spacer RNA at the optimal length for a given Cas13 subtype.
    Extends or trims the default spacer using the junction context in
    fusion_targets.json to match the subtype's preferred spacer length."""
    desired_len = _SUBTYPE_SPACER_LEN.get(subtype, _DEFAULT_SPACER_LEN)
    base_spacer = DUMMY_SPACER_RNA

    if desired_len == len(base_spacer):
        return base_spacer

    info = _TARGET_INFO
    junction_dna = info.get("junction_dna", "")
    if not junction_dna:
        return base_spacer

    junction_rna = junction_dna.upper().replace("T", "U")
    rc_junction = junction_rna[::-1].translate(_RNA_COMPLEMENT)

    if desired_len > len(base_spacer) and len(rc_junction) >= desired_len:
        return rc_junction[:desired_len]
    elif desired_len < len(base_spacer):
        return base_spacer[:desired_len]
    return base_spacer


def get_target_for_spacer(spacer_rna: str) -> str:
    """Build the target RNA (reverse complement + polyA flanks) for a given spacer."""
    return _spacer_to_target(spacer_rna)


def assemble_crrna(dr: str, spacer: str, subtype: str = "unknown") -> str:
    """Assemble the mature crRNA from direct repeat and spacer, respecting
    the subtype-specific orientation:
      - Cas13a/d/X/Y and unknown: 5'-DR-spacer-3'
      - Cas13b:                   5'-spacer-DR-3'
    """
    subtype_lower = (subtype or "unknown").lower()
    if subtype_lower == "cas13b":
        return spacer + dr
    return dr + spacer


def generate_offtarget_sequences(target_rna, num_scrambled=1, num_mismatch=1, seed=None):
    """
    Auto-generate off-target RNA sequences from the target (legacy).
    target_rna: the 24-nt target region (spacer complement).
    Returns list of off-target RNA strings with same flanking AAAAAA structure.
    """
    if seed is not None:
        random.seed(seed)
    offtargets = []
    for _ in range(num_scrambled):
        chars = list(target_rna)
        random.shuffle(chars)
        offtargets.append("AAAAAA" + "".join(chars) + "AAAAAA")
    for _ in range(num_mismatch):
        chars = list(target_rna)
        n_subs = random.randint(2, 3)
        subs = random.sample(range(len(chars)), n_subs)
        for i in subs:
            old = chars[i]
            choices = [c for c in "ACGU" if c != old]
            chars[i] = random.choice(choices)
        offtargets.append("AAAAAA" + "".join(chars) + "AAAAAA")
    return offtargets


def generate_mismatch_sequences(target_rna, mismatch_counts=(1, 2, 3), num_per_count=1, seed=None):
    """
    Generate off-target RNAs with exactly 1, 2, or 3 mismatches.
    Returns list of (rna_string, mismatch_count) tuples.
    Activity at greater mismatch count is worse (less specific) and will be penalized harder.
    """
    if seed is not None:
        random.seed(seed)
    results = []
    for n_mismatch in mismatch_counts:
        for _ in range(num_per_count):
            chars = list(target_rna)
            if n_mismatch > len(chars):
                continue
            subs = random.sample(range(len(chars)), n_mismatch)
            for i in subs:
                old = chars[i]
                choices = [c for c in "ACGU" if c != old]
                chars[i] = random.choice(choices)
            rna = "AAAAAA" + "".join(chars) + "AAAAAA"
            results.append((rna, n_mismatch))
    return results


def generate_offtarget_json(variant_fasta, crrna_lookup_id, metadata_path, off_target_rna, out_dir, suffix=""):
    """
    Generates JSON payload for off-target specificity test (protein + crRNA + off_target_RNA).
    suffix: optional unique suffix for filename (e.g. "0", "1").
    """
    os.makedirs(out_dir, exist_ok=True)
    with open(variant_fasta, 'r') as f:
        protein_seq = "".join([l.strip() for l in f if not l.startswith(">")])
    variant_id = os.path.basename(variant_fasta).replace(".fasta", "")
    with open(metadata_path, 'r') as f:
        metadata = json.load(f)
    baseline_data = metadata.get(crrna_lookup_id)
    if not baseline_data:
        raise ValueError(f"crRNA lookup ID {crrna_lookup_id} not found in metadata.")
    subtype = baseline_data.get("subtype", "unknown")
    spacer = get_spacer_for_subtype(subtype)
    crrna_seq = assemble_crrna(baseline_data["crRNA_repeat_used"], spacer, subtype)
    # Protenix expects proteinChain/rnaSequence (not protein/rna)
    payload = [{
        "name": f"{variant_id}_offtarget_{suffix}",
        "sequences": [
            {"proteinChain": {"sequence": protein_seq, "count": 1}},
            {"rnaSequence": {"sequence": crrna_seq, "count": 1}},
            {"rnaSequence": {"sequence": off_target_rna, "count": 1}}
        ]
    }]
    path = os.path.join(out_dir, f"{variant_id}_offtarget_{suffix}.json")
    with open(path, 'w') as f:
        json.dump(payload, f, indent=2)
    return path


def generate_evaluation_jsons(variant_fasta, baseline_id, metadata_path, out_dir, crrna_lookup_id=None):
    """
    Takes a newly mutated Cas13 FASTA, pairs it with the native crRNA from the baseline,
    and generates the OFF-state and ON-state JSONs required by Protenix.
    crrna_lookup_id: if set, use for metadata lookup instead of baseline_id (for evolved baselines).
    """
    os.makedirs(out_dir, exist_ok=True)
    lookup_id = crrna_lookup_id if crrna_lookup_id is not None else baseline_id

    # 1. Read the mutated protein sequence
    with open(variant_fasta, 'r') as f:
        lines = f.readlines()
        protein_seq = "".join([l.strip() for l in lines if not l.startswith(">")])
        variant_id = os.path.basename(variant_fasta).replace(".fasta", "")

    # 2. Retrieve the native crRNA used for this baseline
    if not os.path.exists(metadata_path):
        raise ValueError(f"Metadata file not found: {metadata_path}")
    with open(metadata_path, 'r') as f:
        metadata = json.load(f)

    baseline_data = metadata.get(lookup_id)
    if not baseline_data:
        raise ValueError(f"crRNA lookup ID {lookup_id} not found in metadata.")

    subtype = baseline_data.get("subtype", "unknown")
    spacer = get_spacer_for_subtype(subtype)
    target_rna = get_target_for_spacer(spacer)
    crrna_seq = assemble_crrna(baseline_data["crRNA_repeat_used"], spacer, subtype)
    
    # 3. Construct OFF State Payload (Dormant - No Target)
    # Protenix expects proteinChain/rnaSequence (not protein/rna)
    off_payload = [{
        "name": f"{variant_id}_OFF",
        "sequences": [
            {"proteinChain": {"sequence": protein_seq, "count": 1}},
            {"rnaSequence": {"sequence": crrna_seq, "count": 1}}
        ]
    }]
    
    # 4. Construct ON State Payload (Triggered - Target Bound)
    on_payload = [{
        "name": f"{variant_id}_ON",
        "sequences": [
            {"proteinChain": {"sequence": protein_seq, "count": 1}},
            {"rnaSequence": {"sequence": crrna_seq, "count": 1}},
            {"rnaSequence": {"sequence": target_rna, "count": 1}}
        ]
    }]
    
    off_json_path = os.path.join(out_dir, f"{variant_id}_OFF.json")
    on_json_path = os.path.join(out_dir, f"{variant_id}_ON.json")
    
    with open(off_json_path, 'w') as f: json.dump(off_payload, f, indent=2)
    with open(on_json_path, 'w') as f: json.dump(on_payload, f, indent=2)
        
    return off_json_path, on_json_path

def _find_cached_outputs(pred_dir):
    """Return (structure_path, summary_path) if both exist under pred_dir, else (None, None)."""
    if not os.path.isdir(pred_dir):
        return None, None
    structure_files = glob.glob(os.path.join(pred_dir, "**/*.cif"), recursive=True)
    if not structure_files:
        structure_files = glob.glob(os.path.join(pred_dir, "**/*.pdb"), recursive=True)
    summary_files = glob.glob(os.path.join(pred_dir, "**/*_summary*.json"), recursive=True)
    if not summary_files:
        summary_files = glob.glob(os.path.join(pred_dir, "*_summary*.json"))
    if structure_files and summary_files:
        return structure_files[0], summary_files[0]
    return None, None


def _find_structure_and_summary(pred_dir):
    """Locate structure (CIF/PDB) and summary JSON under a prediction directory."""
    if not os.path.isdir(pred_dir):
        return None, None
    structure_files = glob.glob(os.path.join(pred_dir, "**/*.cif"), recursive=True)
    if not structure_files:
        structure_files = glob.glob(os.path.join(pred_dir, "*.cif"))
    if not structure_files:
        structure_files = glob.glob(os.path.join(pred_dir, "**/*.pdb"), recursive=True)
    if not structure_files:
        structure_files = glob.glob(os.path.join(pred_dir, "*.pdb"))
    summary_files = glob.glob(os.path.join(pred_dir, "**/*_summary*.json"), recursive=True)
    if not summary_files:
        summary_files = glob.glob(os.path.join(pred_dir, "*_summary*.json"))
    if not summary_files:
        summary_files = glob.glob(os.path.join(pred_dir, "**/*_confidence*.json"), recursive=True)
    if structure_files and summary_files:
        return structure_files[0], summary_files[0]
    return None, None


def _model_name_for_tier(model_tier, engine):
    """Return the model name string for the given tier and engine."""
    if model_tier == "mini":
        if engine == "cattle-prod":
            return "cattle_prod_mini_default_v0.5.0"
        return "protenix_mini_default_v0.5.0"
    if engine == "cattle-prod":
        return os.environ.get("CATTLE_PROD_BASE_MODEL", "cattle_prod_base_default_v1.0.0")
    return os.environ.get("PROTENIX_BASE_MODEL", "protenix_base_default_v1.0.0")


def _resolve_cattle_prod_checkpoint(model_tier):
    """
    Resolve checkpoint directory for cattle-prod and fail loudly if missing.
    Expected env vars:
      - CATTLE_PROD_MINI_CKPT
      - CATTLE_PROD_BASE_CKPT
    """
    env_name = "CATTLE_PROD_MINI_CKPT" if model_tier == "mini" else "CATTLE_PROD_BASE_CKPT"
    ckpt = os.environ.get(env_name, "").strip()
    if not ckpt:
        raise RuntimeError(
            f"{env_name} is required for cattle-prod inference to avoid silent featurize-only mode. "
            f"Set {env_name} to a directory containing model.safetensors."
        )
    if not os.path.isdir(ckpt):
        raise RuntimeError(f"{env_name} path is not a directory: {ckpt}")
    st_path = os.path.join(ckpt, "model.safetensors")
    if not os.path.isfile(st_path):
        raise RuntimeError(f"{env_name} missing model.safetensors: {st_path}")
    return ckpt


def _summary_looks_uninitialized(summary_json_path):
    """
    Detect common silent-featurize signature in confidence summary.
    """
    if not summary_json_path or not os.path.isfile(summary_json_path):
        return False
    try:
        with open(summary_json_path, "r") as f:
            data = json.load(f)
    except Exception:
        return False
    iptm = float(data.get("iptm", 0.0) or 0.0)
    ptm = float(data.get("ptm", 0.0) or 0.0)
    af2_ig = float(data.get("af2_ig", data.get("af2_ig_score", 0.0)) or 0.0)
    ranking = float(data.get("ranking_score", 0.0) or 0.0)
    return iptm == 0.0 and ptm == 0.0 and af2_ig == 0.0 and ranking == 0.0


def _run_pred_command(cmd, engine, base_name, model_tier, base_fallback=None):
    """
    Execute prediction command with engine-specific fallbacks.
    """
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0 and engine == "protenix" and "No such command" in (result.stderr or ""):
        cmd = [
            cmd[0], "predict",
            "--input", cmd[3],  # predict_input
            "--out_dir", cmd[5],  # out_dir
            "--model_name", cmd[7],  # model_name
            "--use_msa", "true",
            "--use_default_params", "true",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0 and model_tier == "base" and base_fallback:
        if "not supported" in (result.stderr or ""):
            log.warning(f"Base model {cmd[7]} not supported; falling back to {base_fallback}")
            cmd = [base_fallback if c == cmd[7] else c for c in cmd]
            result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
    return cmd


def _run_msa_step(json_path, out_dir, base_name, seqres_db_path, engine_bin):
    """Run MSA search step. Returns predict_input_path.

    Always enables MSA for maximum prediction robustness.  The MSA
    enrichment step is attempted (it may add alignment data to the
    JSON); if it fails the original input is used but use_msa remains
    True so the model still activates its MSA pathway.
    """
    msa_dir = os.path.join(out_dir, f"{base_name}_msa")
    os.makedirs(msa_dir, exist_ok=True)
    json_dir = os.path.dirname(os.path.abspath(json_path))
    msa_output_primary = os.path.join(json_dir, f"{base_name}-update-msa.json")
    msa_output_fallback = os.path.join(msa_dir, os.path.basename(json_path))

    for msa_candidate in [msa_output_primary, msa_output_fallback]:
        if os.path.exists(msa_candidate):
            log.info(f"  Reusing cached MSA output for {base_name}")
            return msa_candidate

    try:
        msa_cmd = [engine_bin, "msa", "--input", json_path, "--out_dir", msa_dir]
        if seqres_db_path and os.path.isdir(seqres_db_path):
            msa_cmd.extend(["--db_dir", seqres_db_path])
        subprocess.run(msa_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for msa_candidate in [msa_output_primary, msa_output_fallback]:
            if os.path.exists(msa_candidate):
                return msa_candidate
    except (subprocess.CalledProcessError, FileNotFoundError):
        log.debug(f"  MSA enrichment unavailable for {base_name}; proceeding with use_msa=true on raw input")
    return json_path


def run_protenix_inference(json_path, out_dir, model_tier="mini", seqres_db_path=None):
    """
    Executes the structure prediction engine (cattle-prod or protenix).
    Uses EVAL_CMD env var or auto-detects cattle-prod on PATH.
    model_tier="mini" for Script 4 (Fast Filter)
    model_tier="base" for Script 5 (High Fidelity Oracle)
    seqres_db_path: if set and path exists, runs MSA first for better quality.
    Caches results: skips inference when structure + summary already exist.
    Mini tier skips MSA for speed; base tier runs MSA for accuracy.
    """
    os.makedirs(out_dir, exist_ok=True)
    base_name = os.path.basename(json_path).replace(".json", "")
    tier_label = "mini" if model_tier == "mini" else "base"
    engine_bin, engine = EVAL_BIN, EVAL_ENGINE

    pred_dir = os.path.join(out_dir, base_name)
    cached_struct, cached_summary = _find_cached_outputs(pred_dir)
    if cached_struct and cached_summary:
        log.info(f"Reusing cached {engine} {tier_label} output for {base_name}")
        return cached_struct, cached_summary

    log.info(f"Starting {engine} {tier_label} inference for {base_name} (this may take several minutes)...")
    log.info(f"  eval_decision: engine={engine} strict_mode={_CATTLE_PROD_STRICT}")

    use_msa = model_tier != "mini"
    if use_msa:
        predict_input = _run_msa_step(json_path, out_dir, base_name, seqres_db_path, engine_bin)
    else:
        predict_input = json_path

    model_name = _model_name_for_tier(model_tier, engine)
    checkpoint_dir = None
    if engine == "cattle-prod":
        checkpoint_dir = _resolve_cattle_prod_checkpoint(model_tier)
    base_fallback = None
    if model_tier == "base":
        if engine == "cattle-prod":
            base_fallback = "cattle_prod_base_default_v0.5.0"
        else:
            base_fallback = "protenix_base_default_v0.5.0"

    msa_flag = "true" if use_msa else "false"
    cmd = [
        engine_bin, "pred",
        "-i", predict_input,
        "-o", out_dir,
        "-n", model_name,
        "--use_msa", msa_flag,
        "--use_default_params", "true",
    ]
    if checkpoint_dir:
        cmd.extend(["--checkpoint", checkpoint_dir])

    try:
        _run_pred_command(cmd, engine, base_name, model_tier, base_fallback=base_fallback)
    except subprocess.CalledProcessError as e:
        if engine == "cattle-prod" and not _CATTLE_PROD_STRICT:
            fallback_bin = shutil.which("protenix")
            if fallback_bin:
                log.warning(
                    f"cattle-prod failed for {base_name}, strict mode disabled; "
                    f"falling back to protenix ({fallback_bin})"
                )
                fallback_model = _model_name_for_tier(model_tier, "protenix")
                fb_base_fallback = "protenix_base_default_v0.5.0" if model_tier == "base" else None
                fallback_cmd = [
                    fallback_bin, "pred",
                    "-i", predict_input,
                    "-o", out_dir,
                    "-n", fallback_model,
                    "--use_msa", "true",
                    "--use_default_params", "true",
                ]
                _run_pred_command(
                    fallback_cmd,
                    "protenix",
                    base_name,
                    model_tier,
                    base_fallback=fb_base_fallback,
                )
                engine = "protenix-fallback"
            else:
                print(f"{engine} evaluation failed for {base_name}:\n{e.stderr}")
                raise
        else:
            print(f"{engine} evaluation failed for {base_name}:\n{e.stderr}")
            raise

    struct_path, summary_path = _find_structure_and_summary(pred_dir)
    if not struct_path or not summary_path:
        raise FileNotFoundError(f"{engine} outputs not generated for {base_name}")
    if engine == "cattle-prod" and _summary_looks_uninitialized(summary_path):
        raise RuntimeError(
            f"cattle-prod produced zeroed confidence metrics for {base_name}; "
            f"this usually indicates missing/invalid checkpoint. "
            f"mini_ckpt={os.environ.get('CATTLE_PROD_MINI_CKPT','')} "
            f"base_ckpt={os.environ.get('CATTLE_PROD_BASE_CKPT','')}"
        )

    return struct_path, summary_path