import os
import json
import subprocess
import glob
import random

# RNA Constants (Must match what we defined in 01_parse_and_annotate.py)
DUMMY_SPACER_RNA = "GUCGACUGACGUACGUACGUACGU"
TARGET_REGION = "ACGUACGUACGUACGUCAGUCGAC"  # 24-nt spacer complement
DUMMY_TARGET_RNA = "AAAAAA" + TARGET_REGION + "AAAAAA"


def generate_offtarget_sequences(target_rna, num_scrambled=1, num_mismatch=1, seed=None):
    """
    Auto-generate off-target RNA sequences from the target.
    target_rna: the 24-nt target region (spacer complement).
    Returns list of off-target RNA strings with same flanking AAAAAA structure.
    """
    if seed is not None:
        random.seed(seed)
    offtargets = []
    # Scrambled: random shuffle preserving composition
    for _ in range(num_scrambled):
        chars = list(target_rna)
        random.shuffle(chars)
        offtargets.append("AAAAAA" + "".join(chars) + "AAAAAA")
    # Mismatch: 2-3 random substitutions
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
    crrna_seq = baseline_data["crRNA_repeat_used"] + DUMMY_SPACER_RNA
    payload = [{
        "name": f"{variant_id}_offtarget_{suffix}",
        "sequences": [
            {"protein": {"id": "A", "sequence": protein_seq}},
            {"rna": {"id": "B", "sequence": crrna_seq}},
            {"rna": {"id": "C", "sequence": off_target_rna}}
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
    with open(metadata_path, 'r') as f:
        metadata = json.load(f)

    baseline_data = metadata.get(lookup_id)
    if not baseline_data:
        raise ValueError(f"crRNA lookup ID {lookup_id} not found in metadata.")

    crrna_seq = baseline_data["crRNA_repeat_used"] + DUMMY_SPACER_RNA
    
    # 3. Construct OFF State Payload (Dormant - No Target)
    off_payload = [{
        "name": f"{variant_id}_OFF",
        "sequences": [
            {"protein": {"id": "A", "sequence": protein_seq}},
            {"rna": {"id": "B", "sequence": crrna_seq}}
        ]
    }]
    
    # 4. Construct ON State Payload (Triggered - Target Bound)
    on_payload = [{
        "name": f"{variant_id}_ON",
        "sequences": [
            {"protein": {"id": "A", "sequence": protein_seq}},
            {"rna": {"id": "B", "sequence": crrna_seq}},
            {"rna": {"id": "C", "sequence": DUMMY_TARGET_RNA}}
        ]
    }]
    
    off_json_path = os.path.join(out_dir, f"{variant_id}_OFF.json")
    on_json_path = os.path.join(out_dir, f"{variant_id}_ON.json")
    
    with open(off_json_path, 'w') as f: json.dump(off_payload, f, indent=2)
    with open(on_json_path, 'w') as f: json.dump(on_payload, f, indent=2)
        
    return off_json_path, on_json_path

def run_protenix_inference(json_path, out_dir, model_tier="mini", seqres_db_path=None):
    """
    Executes the Protenix CLI.
    model_tier="mini" for Script 4 (Fast Filter)
    model_tier="base" for Script 5 (High Fidelity Oracle)
    seqres_db_path: if set and path exists, runs inputprep first for better MSA quality (matches Phase 1).
    """
    os.makedirs(out_dir, exist_ok=True)
    base_name = os.path.basename(json_path).replace(".json", "")
    predict_input = json_path

    # Optional: run inputprep first (matches Phase 1 workflow, improves MSA quality)
    if seqres_db_path and os.path.isdir(seqres_db_path):
        prep_dir = os.path.join(out_dir, f"{base_name}_prep")
        os.makedirs(prep_dir, exist_ok=True)
        prep_cmd = [
            "protenix", "inputprep",
            "--input", json_path,
            "--out_dir", prep_dir,
            "--seqres_database_path", seqres_db_path,
        ]
        try:
            subprocess.run(prep_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            # Protenix inputprep typically keeps same filename in output dir
            prep_output = os.path.join(prep_dir, os.path.basename(json_path))
            if os.path.exists(prep_output):
                predict_input = prep_output
        except subprocess.CalledProcessError as e:
            pass  # Fall back to raw JSON if inputprep fails

    if model_tier == "mini":
        model_name = "protenix_mini_default_v0.5.0"
    else:
        model_name = "protenix_base_default_v1.0.0"

    # Construct the RunPod-optimized CLI command
    cmd = [
        "protenix", "predict",
        "--input", predict_input,
        "--out_dir", out_dir,
        "--model_name", model_name,
        "--dtype", "bf16",
        "--enable_cache", "true",
        "--enable_fusion", "true",
        "--trimul_kernel", "true"
    ]
    
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Protenix Evaluation Failed for {base_name}:\n{e.stderr}")
        raise
        
    # Locate outputs
    pred_dir = os.path.join(out_dir, base_name)
    pdb_files = glob.glob(os.path.join(pred_dir, "*.pdb"))
    summary_files = glob.glob(os.path.join(pred_dir, "*_summary.json"))
    
    if not pdb_files or not summary_files:
        raise FileNotFoundError(f"Protenix outputs not generated for {base_name}")
        
    return pdb_files[0], summary_files[0]