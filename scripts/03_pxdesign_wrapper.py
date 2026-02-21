import subprocess
import json
import os
import glob

def generate_frozen_rec_config(metadata_db_path, variant_id, out_dir, metadata_override=None):
    """
    Reads the SQLite/JSON metadata to find the REC lobe boundaries,
    then generates a freeze.json config file for PXDesign.
    metadata_override: optional dict for evolved variants not in metadata file.
    """
    if metadata_override and variant_id in metadata_override:
        variant_data = metadata_override[variant_id]
    else:
        with open(metadata_db_path, 'r') as f:
            metadata = json.load(f)
        variant_data = metadata.get(variant_id)
        if not variant_data:
            raise ValueError(f"Variant {variant_id} not found in metadata.")
        
    # Assume REC lobe starts at 0 and ends right before HEPN1
    hepn1_start = variant_data["domains"]["HEPN1"]["start"]
    rec_end = max(0, hepn1_start - 10) # 10aa buffer
    
    # Create the residue string for PXDesign (e.g., "A1-A350")
    freeze_string = f"A1-A{rec_end}"
    
    freeze_config = {
        "freeze_residues": freeze_string,
        "description": "Strict REC Lobe constraint"
    }
    
    os.makedirs(out_dir, exist_ok=True)
    freeze_path = os.path.join(out_dir, f"{variant_id}_freeze.json")
    with open(freeze_path, 'w') as f:
        json.dump(freeze_config, f, indent=2)
        
    return freeze_path

def run_pxdesign_generation(baseline_pdb, variant_id, metadata_path, bias_json_path, output_dir, variant_count=50, metadata_override=None):
    """
    Executes the physical PXDesign binary on the RunPod cluster, applying 
    both the structural freeze constraints and the Active Learning bias matrix.
    metadata_override: optional dict for evolved variants not in metadata file.
    """
    print(f"[PXDesign] Generating {variant_count} mutations for {variant_id}...")
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Generate the dynamic freeze config to protect the crRNA pocket
    freeze_json_path = generate_frozen_rec_config(metadata_path, variant_id, output_dir, metadata_override)
    
    # 2. Construct the production CLI command
    cmd = [
        "pxdesign",
        "--input_pdb", baseline_pdb,
        "--freeze_json", freeze_json_path,
        "--num_outputs", str(variant_count),
        "--out_dir", output_dir,
        "--model_type", "pxdesign-d" # Utilize the diffusion transformer backbone
    ]
    
    # 3. Apply the Reinforcement Learning Bias Matrix if the Gym provided one
    if bias_json_path and os.path.exists(bias_json_path):
        print(f"  -> Injecting RL Gym Bias Matrix: {os.path.basename(bias_json_path)}")
        cmd.extend(["--bias_aa_json", bias_json_path])
        
    # 4. Execute the command on the GPU
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print(f"PXDesign Execution Failed:\n{e.stderr}")
        raise
        
    # 5. Return the newly generated FASTA files
    generated_fastas = glob.glob(os.path.join(output_dir, f"{variant_id}_variant_*.fasta"))
    return generated_fastas