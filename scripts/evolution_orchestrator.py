import os
import json
import re
import shutil
import importlib.util
from glob import glob
import numpy as np

# --- Import our Production Functional Wrappers ---
# Load pxdesign_wrapper from 03_pxdesign_wrapper.py (module names can't start with digits)
_spec = importlib.util.spec_from_file_location(
    "pxdesign_wrapper",
    os.path.join(os.path.dirname(__file__), "03_pxdesign_wrapper.py"),
)
pxdesign_wrapper = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pxdesign_wrapper)

from utils.protenix_eval import (
    generate_evaluation_jsons,
    run_protenix_inference,
    generate_offtarget_sequences,
    generate_offtarget_json,
    TARGET_REGION,
    DUMMY_SPACER_RNA,
)
from utils.pdb_kinematics import calculate_hepn_shift, extract_protenix_scores

# --- Configuration ---
METADATA_FILE = "../metadata/variant_domain_metadata.json"
BASE_JSON_DIR = "../jsons"
PHASE1_PDB_DIR = "../outputs/phase1_screening"  # Where Script 2 saved the initial PDBs
GENERATION_DIR = "../outputs/generation_queue"
FAST_EVAL_DIR = "../outputs/fast_eval"
HIGH_FIDELITY_DIR = "../outputs/high_fidelity_scoring"
FINAL_HITS_DIR = "../outputs/optimized_switches"
GYM_DIR = "../outputs/rl_gym_data"
# Optional: set to path to databases/ for Protenix inputprep (improves MSA quality)
SEQRES_DB_PATH = os.path.join(os.path.dirname(__file__), "..", "databases")

# --- Biophysical Thresholds ---
# The R(phi)X3H motifs must be far apart in the OFF state, and snap together in the ON state.
MIN_OFF_DISTANCE = 25.0  # Ångströms
MAX_ON_DISTANCE = 12.0   # Ångströms
MIN_IPTM_SCORE = 0.85
MIN_AF2_IG_SCORE = 0.80
# --- Evolution Loop Config ---
MAX_GENERATIONS = 20
NUM_OFFTARGET_TESTS = 2  # scrambled + mismatch variants for specificity
SPECIFICITY_PENALTY_WEIGHT = 0.5
NUM_INITIAL_LINEAGES = 5

def compute_fitness(off_dist, on_dist, iptm_score, af2_ig_score, is_full_ternary=False, offtarget_min_dist=None):
    """
    Composite fitness for ranking variants.
    Higher is better. Incorporates specificity penalty when offtarget_min_dist < MIN_OFF_DISTANCE.
    """
    multiplier = 2.0 if is_full_ternary else 1.0
    shift = off_dist - on_dist
    fitness = (shift - (MIN_OFF_DISTANCE - MAX_ON_DISTANCE)) + ((iptm_score - 0.7) * 50)
    fitness *= multiplier
    if offtarget_min_dist is not None and offtarget_min_dist < MIN_OFF_DISTANCE:
        fitness -= SPECIFICITY_PENALTY_WEIGHT * (MIN_OFF_DISTANCE - offtarget_min_dist)
    return fitness


class EvolutionGym:
    """Active Learning Environment for Directed Evolution."""
    def __init__(self):
        os.makedirs(GYM_DIR, exist_ok=True)
        self.mutation_weights = {}
        self.generation_history = []

    def register_evaluation(self, variant_id, mutations, off_dist, on_dist, iptm_score, is_full_ternary=False, offtarget_min_dist=None):
        """Records the performance of a variant's specific mutations."""
        fitness = compute_fitness(off_dist, on_dist, iptm_score, 0.0, is_full_ternary, offtarget_min_dist)

        self.generation_history.append({
            "variant": variant_id, "fitness": fitness, "mutations": mutations
        })

        for mut in mutations:
            if mut not in self.mutation_weights:
                self.mutation_weights[mut] = 0.0
            self.mutation_weights[mut] = (self.mutation_weights[mut] * 0.5) + (fitness * 0.5)

    def generate_mpnn_bias_matrix(self, generation_num):
        """Converts weights into a physical bias matrix for PXDesign/ProteinMPNN."""
        bias_matrix = {}
        for mut, weight in self.mutation_weights.items():
            parts = mut.split('_')
            if len(parts) >= 2:
                pos = parts[0]
                aa = parts[1] if parts[1] != "del" else None
                if aa and len(aa) == 1:  # Skip 'del', only single-letter AAs for bias
                    if pos not in bias_matrix:
                        bias_matrix[pos] = {}
                    bias_matrix[pos][aa] = float(np.clip(weight / 10.0, -5.0, 5.0))
            
        bias_file = os.path.join(GYM_DIR, f"mpnn_bias_gen_{generation_num}.json")
        with open(bias_file, 'w') as f:
            json.dump(bias_matrix, f, indent=2)
        return bias_file

def build_metadata_override_for_evolved(baseline_id, baseline_fasta_path, crrna_lookup_id, domain_metadata):
    """Build metadata override dict for evolved variants not in variant_domain_metadata.json."""
    with open(baseline_fasta_path, 'r') as f:
        seq = "".join([l.strip() for l in f if not l.startswith(">")])
    motif = re.compile(r'R[AILMFVWY][A-Z]{3}H')
    matches = list(motif.finditer(seq))
    if len(matches) < 2:
        return None
    hepn1_center = matches[0].start()
    hepn2_center = matches[-1].start()
    parent_data = domain_metadata.get(crrna_lookup_id)
    if not parent_data:
        return None
    return {
        baseline_id: {
            "sequence_length": len(seq),
            "domains": {
                "HEPN1": {"start": max(0, hepn1_center - 30), "end": hepn1_center + 80},
                "HEPN2": {"start": max(hepn1_center + 80, hepn2_center - 30), "end": hepn2_center + 80},
            },
            "crRNA_repeat_used": parent_data["crRNA_repeat_used"],
        }
    }


def save_crrna_for_elite(variant_name, crrna_lookup_id, domain_metadata):
    """Saves the crRNA sequence (repeat + spacer) for an elite switch to FINAL_HITS_DIR."""
    parent = domain_metadata.get(crrna_lookup_id)
    if not parent:
        return
    crrna_seq = parent["crRNA_repeat_used"] + DUMMY_SPACER_RNA
    os.makedirs(FINAL_HITS_DIR, exist_ok=True)
    crrna_path = os.path.join(FINAL_HITS_DIR, f"{variant_name}_crRNA.fasta")
    with open(crrna_path, 'w') as f:
        f.write(f">{variant_name}_crRNA\n")
        f.write(f"{crrna_seq}\n")


def get_catalytic_histidine_indices(fasta_path):
    """Parses a FASTA to find the exact 1-based indices of the two catalytic Histidines.
    Uses only the first sequence if the FASTA contains multiple entries."""
    with open(fasta_path, 'r') as f:
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_lines:
                    break  # Stop after first sequence
                continue
            seq_lines.append(line)
        seq = "".join(seq_lines)

    motif = re.compile(r'R[AILMFVWY][A-Z]{3}H')
    matches = list(motif.finditer(seq))
    if len(matches) < 2:
        return None, None
        
    # Matches give the index of 'R'. The 'H' is 5 residues later.
    # Add 1 because Biopython PDB parsing uses 1-based indexing.
    return matches[0].start() + 6, matches[-1].start() + 6

def extract_mutations(baseline_id, variant_fasta, baseline_fasta_path=None):
    """Compares the new variant against the original sequence to map the mutations.
    Handles length mismatches (indels) by padding the shorter sequence.
    baseline_fasta_path: if provided, load baseline sequence from FASTA (for evolved baselines)."""
    if baseline_fasta_path and os.path.exists(baseline_fasta_path):
        with open(baseline_fasta_path, 'r') as f:
            baseline_seq = "".join([l.strip() for l in f if not l.startswith(">")])
    else:
        base_json = os.path.join(BASE_JSON_DIR, f"{baseline_id}.json")
        with open(base_json, 'r') as f:
            data = json.load(f)
        baseline_seq = data[0]["sequences"][0]["protein"]["sequence"]

    with open(variant_fasta, 'r') as f:
        v_seq = "".join([l.strip() for l in f.readlines() if not l.startswith(">")])

    # Pad shorter sequence so we can detect all substitutions and terminal indels
    max_len = max(len(baseline_seq), len(v_seq))
    b_padded = baseline_seq.ljust(max_len, "-")
    v_padded = v_seq.ljust(max_len, "-")

    mutations = []
    for i, (b, v) in enumerate(zip(b_padded, v_padded)):
        if b != v:
            if v == "-":
                mutations.append(f"{i+1}_del")
            elif b == "-":
                mutations.append(f"{i+1}_{v}_ins")
            else:
                mutations.append(f"{i+1}_{v}")
    return mutations

def main_evolution_loop():
    print("Initializing SwitchBlade Active Learning Evolution Loop...")
    os.makedirs(FAST_EVAL_DIR, exist_ok=True)
    os.makedirs(HIGH_FIDELITY_DIR, exist_ok=True)

    gym = EvolutionGym()

    with open(METADATA_FILE, 'r') as f:
        domain_metadata = json.load(f)

    # Baseline object: (baseline_id, baseline_pdb_path, baseline_fasta_path, crrna_lookup_id)
    lineage_queue = [
        (bid, None, None, bid)  # Phase 1: fasta_path=None, crrna_lookup_id=baseline_id
        for bid in list(domain_metadata.keys())[:NUM_INITIAL_LINEAGES]
    ]
    if not lineage_queue:
        print("No baselines in metadata. Exiting.")
        return

    baseline = lineage_queue.pop(0)
    baseline_id, baseline_pdb_path, baseline_fasta_path, crrna_lookup_id = baseline

    # Resolve Phase 1 PDB path
    pdbs = glob(f"{PHASE1_PDB_DIR}/{baseline_id}_pred/*.pdb")
    if not pdbs:
        print(f"Warning: No Phase 1 PDB found for {baseline_id}. Trying next lineage...")
        while lineage_queue:
            baseline = lineage_queue.pop(0)
            baseline_id, _, _, crrna_lookup_id = baseline
            pdbs = glob(f"{PHASE1_PDB_DIR}/{baseline_id}_pred/*.pdb")
            if pdbs:
                baseline_pdb_path = pdbs[0]
                baseline = (baseline_id, baseline_pdb_path, None, crrna_lookup_id)
                break
        else:
            print("No valid Phase 1 PDBs found. Exiting.")
            return
    else:
        baseline_pdb_path = pdbs[0]
        baseline = (baseline_id, baseline_pdb_path, baseline_fasta_path, crrna_lookup_id)

    generation_counter = 0
    bias_file = None
    offtarget_seqs = generate_offtarget_sequences(TARGET_REGION, num_scrambled=1, num_mismatch=1, seed=42)

    for generation_counter in range(1, MAX_GENERATIONS + 1):
        baseline_id, baseline_pdb_path, baseline_fasta_path, crrna_lookup_id = baseline

        print(f"\n=======================================================")
        print(f"=== Evolution Generation {generation_counter} | Baseline: {baseline_id} ===")
        print(f"=======================================================")

        metadata_override = None
        if baseline_fasta_path:
            metadata_override = build_metadata_override_for_evolved(baseline_id, baseline_fasta_path, crrna_lookup_id, domain_metadata)
            if not metadata_override:
                print(f"Warning: Could not build metadata for evolved baseline {baseline_id}. Skipping.")
                if lineage_queue:
                    baseline = lineage_queue.pop(0)
                    baseline_id, _, _, crrna_lookup_id = baseline
                    pdbs = glob(f"{PHASE1_PDB_DIR}/{baseline_id}_pred/*.pdb")
                    if pdbs:
                        baseline = (baseline_id, pdbs[0], None, crrna_lookup_id)
                    continue

        # --- SCRIPT 3: Generate Variants (Guided by the Gym) ---
        try:
            new_variants_fastas = pxdesign_wrapper.run_pxdesign_generation(
                baseline_pdb=baseline_pdb_path,
                variant_id=baseline_id,
                metadata_path=METADATA_FILE,
                bias_json_path=bias_file,
                output_dir=os.path.join(GENERATION_DIR, f"gen_{generation_counter}"),
                variant_count=25,
                metadata_override=metadata_override,
            )
        except Exception as e:
            print(f"PXDesign failed: {e}. Trying next lineage...")
            if lineage_queue:
                baseline = lineage_queue.pop(0)
                baseline_id, _, _, crrna_lookup_id = baseline
                pdbs = glob(f"{PHASE1_PDB_DIR}/{baseline_id}_pred/*.pdb")
                if pdbs:
                    baseline = (baseline_id, pdbs[0], None, crrna_lookup_id)
            continue

        if not new_variants_fastas:
            print("No variants generated. Trying next lineage...")
            if lineage_queue:
                baseline = lineage_queue.pop(0)
                baseline_id, _, _, crrna_lookup_id = baseline
                pdbs = glob(f"{PHASE1_PDB_DIR}/{baseline_id}_pred/*.pdb")
                if pdbs:
                    baseline = (baseline_id, pdbs[0], None, crrna_lookup_id)
            continue

        results = []

        for variant_fasta in new_variants_fastas:
            mutations_made = extract_mutations(baseline_id, variant_fasta, baseline_fasta_path)
            variant_name = os.path.basename(variant_fasta).replace(".fasta", "")

            h1_idx, h2_idx = get_catalytic_histidine_indices(variant_fasta)
            if not h1_idx:
                continue

            print(f"  -> Rapid Evaluating {variant_name}...")
            off_json, on_json = generate_evaluation_jsons(
                variant_fasta, baseline_id, METADATA_FILE, FAST_EVAL_DIR, crrna_lookup_id=crrna_lookup_id
            )

            try:
                off_pdb, _ = run_protenix_inference(
                    off_json, FAST_EVAL_DIR, model_tier="mini", seqres_db_path=SEQRES_DB_PATH
                )
                on_pdb, _ = run_protenix_inference(
                    on_json, FAST_EVAL_DIR, model_tier="mini", seqres_db_path=SEQRES_DB_PATH
                )
            except Exception as e:
                print(f"     Protenix failed for {variant_name}: {e}")
                fitness = compute_fitness(0, 999, 0.4, 0, False, None)
                gym.register_evaluation(variant_name, mutations_made, 0, 999, 0.4, False, None)
                results.append((variant_name, variant_fasta, fitness, 0, 999, 0.4, 0.0, None, None))
                continue

            off_dist = calculate_hepn_shift(off_pdb, h1_idx, h2_idx)
            on_dist = calculate_hepn_shift(on_pdb, h1_idx, h2_idx)
            print(f"     [Filter] OFF: {off_dist:.1f}A | ON: {on_dist:.1f}A")
            has_potential = (off_dist >= MIN_OFF_DISTANCE) and (on_dist <= MAX_ON_DISTANCE)

            offtarget_min_dist = None
            hf_pdb_path = None
            iptm, af2_ig = 0.4, 0.0
            true_on_dist = on_dist

            if has_potential:
                print(f"  Filter Passed. Triggering High-Fidelity Ternary Oracle...")

                hf_pdb, hf_summary = run_protenix_inference(
                    on_json, HIGH_FIDELITY_DIR, model_tier="base", seqres_db_path=SEQRES_DB_PATH
                )
                true_on_dist = calculate_hepn_shift(hf_pdb, h1_idx, h2_idx)
                scores = extract_protenix_scores(hf_summary)
                iptm, af2_ig = scores["iptm"], scores["af2_ig"]
                hf_pdb_path = hf_pdb

                # Specificity: off-target tests
                offtarget_dists = []
                for i, ot_rna in enumerate(offtarget_seqs[:NUM_OFFTARGET_TESTS]):
                    try:
                        ot_json = generate_offtarget_json(
                            variant_fasta, crrna_lookup_id, METADATA_FILE, ot_rna, FAST_EVAL_DIR, suffix=str(i)
                        )
                        ot_pdb, _ = run_protenix_inference(
                            ot_json, FAST_EVAL_DIR, model_tier="mini", seqres_db_path=SEQRES_DB_PATH
                        )
                        ot_dist = calculate_hepn_shift(ot_pdb, h1_idx, h2_idx)
                        offtarget_dists.append(ot_dist)
                    except Exception:
                        offtarget_dists.append(MIN_OFF_DISTANCE)  # Assume specific on failure
                offtarget_min_dist = min(offtarget_dists) if offtarget_dists else None

            fitness = compute_fitness(off_dist, true_on_dist, iptm, af2_ig, has_potential, offtarget_min_dist)
            gym.register_evaluation(
                variant_name, mutations_made, off_dist, true_on_dist, iptm, has_potential, offtarget_min_dist
            )
            results.append((variant_name, variant_fasta, fitness, off_dist, true_on_dist, iptm, af2_ig, hf_pdb_path, offtarget_min_dist))

        if not results:
            print("No valid results this generation. Trying next lineage...")
            if lineage_queue:
                baseline = lineage_queue.pop(0)
                baseline_id, _, _, crrna_lookup_id = baseline
                pdbs = glob(f"{PHASE1_PDB_DIR}/{baseline_id}_pred/*.pdb")
                if pdbs:
                    baseline = (baseline_id, pdbs[0], None, crrna_lookup_id)
            continue

        best = max(results, key=lambda r: r[2])
        best_name, best_fasta, best_fitness, best_off, best_on, best_iptm, best_af2_ig, best_hf_pdb, _ = best

        print(f"\n  Best variant: {best_name} (fitness={best_fitness:.1f})")

        if best_iptm >= MIN_IPTM_SCORE and best_af2_ig >= MIN_AF2_IG_SCORE and best_on <= MAX_ON_DISTANCE:
            print(f"  ELITE TERNARY SWITCH FOUND! ipTM: {best_iptm:.3f} | AF2-IG: {best_af2_ig:.3f} | ON-Dist: {best_on:.1f}A")
            os.makedirs(FINAL_HITS_DIR, exist_ok=True)
            shutil.copy(best_fasta, os.path.join(FINAL_HITS_DIR, f"{best_name}_optimal.fasta"))
            if best_hf_pdb:
                shutil.copy(best_hf_pdb, os.path.join(FINAL_HITS_DIR, f"{best_name}_ternary_complex.pdb"))
            save_crrna_for_elite(best_name, crrna_lookup_id, domain_metadata)

        # Use best variant as next baseline (even if not elite)
        os.makedirs(FINAL_HITS_DIR, exist_ok=True)
        best_fasta_dest = os.path.join(FINAL_HITS_DIR, f"{best_name}_optimal.fasta")
        best_pdb_dest = os.path.join(FINAL_HITS_DIR, f"{best_name}_ternary_complex.pdb")

        if best_hf_pdb:
            shutil.copy(best_fasta, best_fasta_dest)
            shutil.copy(best_hf_pdb, best_pdb_dest)
            save_crrna_for_elite(best_name, crrna_lookup_id, domain_metadata)
            baseline = (best_name, best_pdb_dest, best_fasta_dest, crrna_lookup_id)
        elif os.path.exists(best_pdb_dest):
            shutil.copy(best_fasta, best_fasta_dest)
            save_crrna_for_elite(best_name, crrna_lookup_id, domain_metadata)
            baseline = (best_name, best_pdb_dest, best_fasta_dest, crrna_lookup_id)
        else:
            on_json = os.path.join(FAST_EVAL_DIR, f"{best_name}_ON.json")
            if os.path.exists(on_json):
                try:
                    hf_pdb, _ = run_protenix_inference(
                        on_json, HIGH_FIDELITY_DIR, model_tier="base", seqres_db_path=SEQRES_DB_PATH
                    )
                    shutil.copy(best_fasta, best_fasta_dest)
                    shutil.copy(hf_pdb, best_pdb_dest)
                    save_crrna_for_elite(best_name, crrna_lookup_id, domain_metadata)
                    baseline = (best_name, best_pdb_dest, best_fasta_dest, crrna_lookup_id)
                except Exception as e:
                    print(f"  Could not get HF PDB for best variant: {e}. Reusing current baseline.")
                    baseline = (baseline_id, baseline_pdb_path, baseline_fasta_path, crrna_lookup_id)
            else:
                baseline = (baseline_id, baseline_pdb_path, baseline_fasta_path, crrna_lookup_id)

        if gym.mutation_weights:
            bias_file = gym.generate_mpnn_bias_matrix(generation_counter)

    print("\nEvolution loop complete.")

if __name__ == "__main__":
    main_evolution_loop()