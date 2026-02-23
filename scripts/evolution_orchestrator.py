import os
import json
import re
import shutil
import gc
import time
import logging
import datetime
import importlib.util
import numpy as np

# --- Logging ---
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

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
    generate_mismatch_sequences,
    generate_offtarget_json,
    TARGET_REGION,
    DUMMY_SPACER_RNA,
)
from utils.pdb_kinematics import calculate_hepn_shift, extract_protenix_scores, find_structure_files

# --- Configuration ---
RUN_ID = os.environ.get("CASCADE_RUN_ID") or datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
OUTPUTS_RUN = os.path.join(os.path.dirname(__file__), "..", "outputs", f"run_{RUN_ID}")
METADATA_FILE = "../metadata/variant_domain_metadata.json"
BASE_JSON_DIR = "../jsons"
PHASE1_PDB_DIR = "../outputs/phase1_screening"  # Where Script 2 saved the initial PDBs (shared across runs)
GENERATION_DIR = os.path.join(OUTPUTS_RUN, "generation_queue")
FAST_EVAL_DIR = os.path.join(OUTPUTS_RUN, "fast_eval")
HIGH_FIDELITY_DIR = os.path.join(OUTPUTS_RUN, "high_fidelity_scoring")
FINAL_HITS_DIR = os.path.join(OUTPUTS_RUN, "optimized_switches")
GYM_DIR = os.path.join(OUTPUTS_RUN, "rl_gym_data")
RL_TRAINING_DATASET = os.path.join(GYM_DIR, "rl_training_dataset.jsonl")
VALIDATED_IDS_FILE = "../outputs/validated_baseline_ids.txt"  # From validate_crispr_repeats.py; restricts lineage to validated repeats
# Optional: set to path to databases/ for Protenix inputprep (improves MSA quality)
SEQRES_DB_PATH = os.path.join(os.path.dirname(__file__), "..", "databases")


def _get_next_baseline_from_queue(lineage_queue):
    """Pop from queue until we find one with a valid Phase 1 structure. Returns (baseline, lineage_queue) or (None, lineage_queue) if none found."""
    while lineage_queue:
        baseline = lineage_queue.pop(0)
        baseline_id, _, _, crrna_lookup_id = baseline
        structures = find_structure_files(os.path.join(PHASE1_PDB_DIR, f"{baseline_id}_pred"))
        if structures:
            return (baseline_id, structures[0], None, crrna_lookup_id), lineage_queue
        log.warning(f"No Phase 1 structure for {baseline_id}. Skipping to next lineage.")
    return None, lineage_queue


def _load_validated_baseline_ids():
    """Load baseline IDs that passed CRISPR repeat validation. Returns None if file missing (use all)."""
    path = os.path.join(os.path.dirname(__file__), VALIDATED_IDS_FILE)
    if not os.path.isfile(path):
        return None
    ids = set()
    with open(path) as f:
        for line in f:
            bid = line.strip()
            if bid:
                ids.add(bid)
    return ids

# --- Biophysical Thresholds ---
# The R(phi)X3H motifs must be far apart in the OFF state, and snap together in the ON state.
MIN_OFF_DISTANCE = 25.0  # Ångströms
MAX_ON_DISTANCE = 12.0   # Ångströms
MIN_IPTM_SCORE = 0.85
MIN_AF2_IG_SCORE = 0.80
# --- Evolution Loop Config ---
MAX_GENERATIONS = 20
MISMATCH_COUNTS = (1, 2, 3)  # Test 1-, 2-, 3-mismatch off-targets; activity at higher count penalized harder
SPECIFICITY_PENALTY_BASE = 0.3  # Base penalty; scaled by mismatch count (3mm > 2mm > 1mm)
FALLBACK_FITNESS_PENALTY = 5.0  # Penalty when PXDesign stitching fails and baseline is used as fallback
# --- Memory / OOM mitigation (seconds; set to 0 to disable) ---
SLEEP_AFTER_PXDESIGN = 2.0       # Allow CUDA driver to reclaim GPU memory after diffusion
SLEEP_AFTER_PROTENIX_BASE = 2.0  # Allow reclaim after heavy base-model ternary prediction
SLEEP_AFTER_PROTENIX_MINI = 0.0  # Optional: use 0.5 if OOM on back-to-back mini runs

def compute_fitness(off_dist, on_dist, iptm_score, af2_ig_score, is_full_ternary=False, offtarget_by_mismatch=None):
    """
    Composite fitness for ranking variants.
    Primary metric: non-active when unbound (high OFF dist), highly active when bound (low ON dist, high ipTM).
    Specificity: penalize activity on 1-, 2-, 3-mismatch off-targets; greater mismatch activation = harder penalty.
    """
    multiplier = 2.0 if is_full_ternary else 1.0
    shift = off_dist - on_dist
    fitness = (shift - (MIN_OFF_DISTANCE - MAX_ON_DISTANCE)) + ((iptm_score - 0.7) * 50)
    fitness *= multiplier

    # Progressive specificity penalty: activity at 3mm > 2mm > 1mm penalized harder
    if offtarget_by_mismatch:
        for n_mismatch, min_dist in offtarget_by_mismatch.items():
            if min_dist is not None and min_dist < MIN_OFF_DISTANCE:
                penalty = SPECIFICITY_PENALTY_BASE * n_mismatch * (MIN_OFF_DISTANCE - min_dist)
                fitness -= penalty
    return fitness


class EvolutionGym:
    """Active Learning Environment for Directed Evolution."""
    def __init__(self):
        os.makedirs(GYM_DIR, exist_ok=True)
        self.mutation_weights = {}
        self.generation_history = []

    def register_evaluation(self, variant_id, mutations, off_dist, on_dist, iptm_score, is_full_ternary=False, offtarget_by_mismatch=None):
        """Records the performance of a variant's specific mutations."""
        fitness = compute_fitness(off_dist, on_dist, iptm_score, 0.0, is_full_ternary, offtarget_by_mismatch)

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


def _read_sequence_from_fasta(fasta_path):
    """Read first sequence from FASTA file."""
    if not fasta_path or not os.path.exists(fasta_path):
        return ""
    with open(fasta_path) as f:
        return "".join(l.strip() for l in f if not l.startswith(">"))


def _read_baseline_sequence(baseline_id, baseline_fasta_path):
    """Get baseline sequence from FASTA or base JSON."""
    if baseline_fasta_path and os.path.exists(baseline_fasta_path):
        return _read_sequence_from_fasta(baseline_fasta_path)
    base_json = os.path.join(BASE_JSON_DIR, f"{baseline_id}.json")
    if os.path.exists(base_json):
        with open(base_json) as f:
            data = json.load(f)
        ent = data[0]["sequences"][0]
        prot = ent.get("proteinChain", ent.get("protein", {}))
        return prot.get("sequence", "")
    return ""


def save_rl_training_record(
    variant_id, variant_fasta, baseline_id, baseline_fasta_path, crrna_lookup_id,
    generation, mutations, fitness, off_dist, on_dist, iptm, af2_ig,
    structure_path, offtarget_by_mismatch, is_elite,
):
    """
    Append a single variant evaluation to rl_training_dataset.jsonl.
    Format is designed for DRAKES/ProteinMPNN post-training: (structure, sequence, reward).
    After post-training MPNN on this data, the fine-tuned model can be plugged back
    into the CASCADE pipeline to bias future designs toward better Cas13 switches.
    """
    os.makedirs(GYM_DIR, exist_ok=True)
    seq = _read_sequence_from_fasta(variant_fasta)
    baseline_seq = _read_baseline_sequence(baseline_id, baseline_fasta_path)
    record = {
        "variant_id": variant_id,
        "generation": generation,
        "baseline_id": baseline_id,
        "crrna_lookup_id": crrna_lookup_id,
        "sequence": seq,
        "baseline_sequence": baseline_seq,
        "mutations": mutations,
        "fitness": float(fitness),
        "off_dist_A": float(off_dist),
        "on_dist_A": float(on_dist),
        "iptm": float(iptm),
        "af2_ig": float(af2_ig),
        "structure_path": structure_path,
        "offtarget_by_mismatch": offtarget_by_mismatch or {},
        "is_elite": bool(is_elite),
    }
    with open(RL_TRAINING_DATASET, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")


def build_metadata_override_for_evolved(baseline_id, baseline_fasta_path, crrna_lookup_id, domain_metadata):
    """Build metadata override dict for evolved variants not in variant_domain_metadata.json."""
    with open(baseline_fasta_path, 'r') as f:
        seq = "".join([l.strip() for l in f if not l.startswith(">")])
    motif = re.compile(r'R.{3,6}H')  # includes Cas13a (REFYH)
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

    motif = re.compile(r'R.{3,6}H')  # includes Cas13a (REFYH)
    matches = list(motif.finditer(seq))
    if len(matches) < 2:
        return None, None
        
    # H is at the end of each match. match.end() gives 1-based index (Biopython PDB uses 1-based).
    return matches[0].end(), matches[-1].end()

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
        baseline_seq = data[0]["sequences"][0].get("proteinChain", data[0]["sequences"][0].get("protein", {}))["sequence"]

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
    log.info("Initializing SwitchBlade Active Learning Evolution Loop...")
    os.makedirs(FAST_EVAL_DIR, exist_ok=True)
    os.makedirs(HIGH_FIDELITY_DIR, exist_ok=True)
    os.makedirs(GYM_DIR, exist_ok=True)
    log.info(f"RL training data will be appended to {RL_TRAINING_DATASET} (see RL_TRAINING_FORMAT.md)")

    gym = EvolutionGym()

    with open(METADATA_FILE, 'r') as f:
        domain_metadata = json.load(f)

    baseline_ids = list(domain_metadata.keys())
    validated = _load_validated_baseline_ids()
    if validated is not None:
        baseline_ids = [b for b in baseline_ids if b in validated]
        log.info(f"Restricting to {len(baseline_ids)} baselines with validated CRISPR repeats (from {VALIDATED_IDS_FILE})")
    else:
        log.info(f"Using all {len(baseline_ids)} baselines (no {VALIDATED_IDS_FILE})")

    # Baseline object: (baseline_id, baseline_pdb_path, baseline_fasta_path, crrna_lookup_id)
    lineage_queue = [
        (bid, None, None, bid)  # Phase 1: fasta_path=None, crrna_lookup_id=baseline_id
        for bid in baseline_ids
    ]
    log.info(f"Queue contains {len(lineage_queue)} lineages ({MAX_GENERATIONS} generations each)")
    if not lineage_queue:
        log.warning("No baselines in metadata. Exiting.")
        return

    baseline, lineage_queue = _get_next_baseline_from_queue(lineage_queue)
    if baseline is None:
        log.warning("No valid Phase 1 structures found for any baseline. Exiting.")
        return

    generation_counter = 0
    bias_file = None
    mismatch_seqs = generate_mismatch_sequences(TARGET_REGION, mismatch_counts=MISMATCH_COUNTS, num_per_count=1, seed=42)

    while True:
        champion = None  # Elitist: highest fitness so far in this lineage; only update when we find better
        for generation_counter in range(1, MAX_GENERATIONS + 1):
            baseline_id, baseline_pdb_path, baseline_fasta_path, crrna_lookup_id = baseline

            log.info("=" * 60)
            log.info(f"Evolution Generation {generation_counter} | Baseline: {baseline_id}")
            log.info("=" * 60)

            metadata_override = None
            if baseline_fasta_path:
                metadata_override = build_metadata_override_for_evolved(baseline_id, baseline_fasta_path, crrna_lookup_id, domain_metadata)
                if not metadata_override:
                    log.warning(f"Could not build metadata for evolved baseline {baseline_id}. Skipping.")
                    if lineage_queue:
                        baseline, lineage_queue = _get_next_baseline_from_queue(lineage_queue)
                        if baseline is not None:
                            continue
                    break

            # --- SCRIPT 3: Generate Variants (Guided by the Gym) ---
            try:
                new_variants_fastas = pxdesign_wrapper.run_pxdesign_generation(
                    baseline_structure=baseline_pdb_path,
                    variant_id=baseline_id,
                    metadata_path=METADATA_FILE,
                    bias_json_path=bias_file,
                    output_dir=os.path.join(GENERATION_DIR, f"gen_{generation_counter}"),
                    variant_count=2,
                    metadata_override=metadata_override,
                    baseline_fasta_path=baseline_fasta_path,
                    base_json_dir=BASE_JSON_DIR,
                )
            except Exception as e:
                log.error(f"PXDesign failed: {e}. Trying next lineage...")
                if lineage_queue:
                    baseline, lineage_queue = _get_next_baseline_from_queue(lineage_queue)
                    if baseline is not None:
                        continue
                break

            if SLEEP_AFTER_PXDESIGN > 0:
                time.sleep(SLEEP_AFTER_PXDESIGN)
            gc.collect()

            if not new_variants_fastas:
                print("No variants generated. Trying next lineage...")
                if lineage_queue:
                    baseline, lineage_queue = _get_next_baseline_from_queue(lineage_queue)
                    if baseline is not None:
                        continue
                break

            results = []

            for variant_fasta in new_variants_fastas:
                mutations_made = extract_mutations(baseline_id, variant_fasta, baseline_fasta_path)
                variant_name = os.path.basename(variant_fasta).replace(".fasta", "")

                h1_idx, h2_idx = get_catalytic_histidine_indices(variant_fasta)
                if not h1_idx:
                    continue

                log.info(f"Evaluating variant {variant_name} (Protenix mini OFF/ON - may take 2-5 min each)...")
                off_json, on_json = generate_evaluation_jsons(
                    variant_fasta, baseline_id, METADATA_FILE, FAST_EVAL_DIR, crrna_lookup_id=crrna_lookup_id
                )

                try:
                    off_pdb, _ = run_protenix_inference(
                        off_json, FAST_EVAL_DIR, model_tier="mini", seqres_db_path=SEQRES_DB_PATH
                    )
                    if SLEEP_AFTER_PROTENIX_MINI > 0:
                        time.sleep(SLEEP_AFTER_PROTENIX_MINI)
                    on_pdb, _ = run_protenix_inference(
                        on_json, FAST_EVAL_DIR, model_tier="mini", seqres_db_path=SEQRES_DB_PATH
                    )
                    if SLEEP_AFTER_PROTENIX_MINI > 0:
                        time.sleep(SLEEP_AFTER_PROTENIX_MINI)
                except Exception as e:
                    log.warning(f"Protenix failed for {variant_name}: {e}")
                    fitness = compute_fitness(0, 999, 0.4, 0, False, None)
                    gym.register_evaluation(variant_name, mutations_made, 0, 999, 0.4, False, None)
                    save_rl_training_record(
                        variant_name, variant_fasta, baseline_id, baseline_fasta_path, crrna_lookup_id,
                        generation_counter, mutations_made, fitness, 0, 999, 0.4, 0.0,
                        None, None, False,
                    )
                    results.append((variant_name, variant_fasta, fitness, 0, 999, 0.4, 0.0, None, None))
                    continue

                off_dist = calculate_hepn_shift(off_pdb, h1_idx, h2_idx)
                on_dist = calculate_hepn_shift(on_pdb, h1_idx, h2_idx)
                print(f"     [Filter] OFF: {off_dist:.1f}A | ON: {on_dist:.1f}A")
                has_potential = (off_dist >= MIN_OFF_DISTANCE) and (on_dist <= MAX_ON_DISTANCE)

                offtarget_by_mismatch = {}
                hf_pdb_path = None
                iptm, af2_ig = 0.4, 0.0
                true_on_dist = on_dist

                if has_potential:
                    log.info("Filter passed. Running Protenix base ternary (may take 10-30 min)...")

                    hf_pdb, hf_summary = run_protenix_inference(
                        on_json, HIGH_FIDELITY_DIR, model_tier="base", seqres_db_path=SEQRES_DB_PATH
                    )
                    if SLEEP_AFTER_PROTENIX_BASE > 0:
                        time.sleep(SLEEP_AFTER_PROTENIX_BASE)
                    gc.collect()
                    true_on_dist = calculate_hepn_shift(hf_pdb, h1_idx, h2_idx)
                    scores = extract_protenix_scores(hf_summary)
                    iptm, af2_ig = scores["iptm"], scores["af2_ig"]
                    hf_pdb_path = hf_pdb

                    # Specificity: 1-, 2-, 3-mismatch off-target tests; activity at higher count penalized harder
                    for i, (ot_rna, n_mismatch) in enumerate(mismatch_seqs):
                        try:
                            ot_json = generate_offtarget_json(
                                variant_fasta, crrna_lookup_id, METADATA_FILE, ot_rna, FAST_EVAL_DIR,
                                suffix=f"{n_mismatch}mm_{i}"
                            )
                            ot_pdb, _ = run_protenix_inference(
                                ot_json, FAST_EVAL_DIR, model_tier="mini", seqres_db_path=SEQRES_DB_PATH
                            )
                            if SLEEP_AFTER_PROTENIX_MINI > 0:
                                time.sleep(SLEEP_AFTER_PROTENIX_MINI)
                            ot_dist = calculate_hepn_shift(ot_pdb, h1_idx, h2_idx)
                            if n_mismatch not in offtarget_by_mismatch:
                                offtarget_by_mismatch[n_mismatch] = ot_dist
                            else:
                                offtarget_by_mismatch[n_mismatch] = min(offtarget_by_mismatch[n_mismatch], ot_dist)
                        except Exception:
                            offtarget_by_mismatch[n_mismatch] = MIN_OFF_DISTANCE  # Assume specific on failure
                    if offtarget_by_mismatch:
                        mm_str = " | ".join(f"{k}mm:{v:.1f}A" for k, v in sorted(offtarget_by_mismatch.items()))
                        print(f"     [Specificity] {mm_str}")

                fitness = compute_fitness(off_dist, true_on_dist, iptm, af2_ig, has_potential, offtarget_by_mismatch or None)
                if "fallback" in variant_name:
                    fitness -= FALLBACK_FITNESS_PENALTY
                    log.info(f"Fallback variant {variant_name}: applying penalty ({FALLBACK_FITNESS_PENALTY})")
                gym.register_evaluation(
                    variant_name, mutations_made, off_dist, true_on_dist, iptm, has_potential, offtarget_by_mismatch or None
                )
                struct_path = hf_pdb_path if hf_pdb_path else on_pdb
                is_elite = (
                    off_dist >= MIN_OFF_DISTANCE and true_on_dist <= MAX_ON_DISTANCE
                    and iptm >= MIN_IPTM_SCORE and af2_ig >= MIN_AF2_IG_SCORE
                )
                save_rl_training_record(
                    variant_name, variant_fasta, baseline_id, baseline_fasta_path, crrna_lookup_id,
                    generation_counter, mutations_made, fitness, off_dist, true_on_dist, iptm, af2_ig,
                    struct_path, offtarget_by_mismatch or None, is_elite,
                )
                results.append((variant_name, variant_fasta, fitness, off_dist, true_on_dist, iptm, af2_ig, hf_pdb_path, offtarget_by_mismatch))
                gc.collect()

            if not results:
                log.warning("No valid results this generation. Trying next lineage...")
                if lineage_queue:
                    baseline, lineage_queue = _get_next_baseline_from_queue(lineage_queue)
                    if baseline is not None:
                        continue
                break

            best = max(results, key=lambda r: r[2])
            best_name, best_fasta, best_fitness, best_off, best_on, best_iptm, best_af2_ig, best_hf_pdb, _ = best

            log.info(f"Best variant this gen: {best_name} (fitness={best_fitness:.1f})")

            os.makedirs(FINAL_HITS_DIR, exist_ok=True)

            # Save all variants that meet the elite threshold to FINAL_HITS_DIR (even if not champion)
            def _is_elite(off, on, iptm, af2_ig):
                return (
                    off >= MIN_OFF_DISTANCE and on <= MAX_ON_DISTANCE
                    and iptm >= MIN_IPTM_SCORE and af2_ig >= MIN_AF2_IG_SCORE
                )

            def _find_existing_structure(final_dir, name):
                base = os.path.join(final_dir, f"{name}_ternary_complex")
                for ext in (".cif", ".pdb"):
                    p = base + ext
                    if os.path.exists(p):
                        return p
                return None

            for r in results:
                v_name, v_fasta, _, v_off, v_on, v_iptm, v_af2_ig, v_hf_pdb, _ = r
                if _is_elite(v_off, v_on, v_iptm, v_af2_ig):
                    log.info(f"ELITE TERNARY SWITCH: {v_name} (ipTM: {v_iptm:.3f} | AF2-IG: {v_af2_ig:.3f} | OFF: {v_off:.1f}A | ON: {v_on:.1f}A)")
                    fasta_dest = os.path.join(FINAL_HITS_DIR, f"{v_name}_optimal.fasta")
                    shutil.copy(v_fasta, fasta_dest)
                    if v_hf_pdb:
                        ext = os.path.splitext(v_hf_pdb)[1]
                        struct_dest = os.path.join(FINAL_HITS_DIR, f"{v_name}_ternary_complex{ext}")
                        shutil.copy(v_hf_pdb, struct_dest)
                    save_crrna_for_elite(v_name, crrna_lookup_id, domain_metadata)

            # Resolve best variant's structure for champion candidate
            best_resolved = False
            best_structure_ext = os.path.splitext(best_hf_pdb)[1] if best_hf_pdb else ".pdb"
            best_structure_dest = os.path.join(FINAL_HITS_DIR, f"{best_name}_ternary_complex{best_structure_ext}")
            best_fasta_dest = os.path.join(FINAL_HITS_DIR, f"{best_name}_optimal.fasta")

            if best_hf_pdb:
                shutil.copy(best_fasta, best_fasta_dest)
                shutil.copy(best_hf_pdb, best_structure_dest)
                best_resolved = True
            elif (existing_structure := _find_existing_structure(FINAL_HITS_DIR, best_name)):
                shutil.copy(best_fasta, best_fasta_dest)
                best_structure_dest = existing_structure
                best_resolved = True
            else:
                on_json = os.path.join(FAST_EVAL_DIR, f"{best_name}_ON.json")
                if os.path.exists(on_json):
                    try:
                        hf_structure, _ = run_protenix_inference(
                            on_json, HIGH_FIDELITY_DIR, model_tier="base", seqres_db_path=SEQRES_DB_PATH
                        )
                        shutil.copy(best_fasta, best_fasta_dest)
                        ext = os.path.splitext(hf_structure)[1] or ".pdb"
                        best_structure_dest = os.path.join(FINAL_HITS_DIR, f"{best_name}_ternary_complex{ext}")
                        shutil.copy(hf_structure, best_structure_dest)
                        best_resolved = True
                    except Exception as e:
                        print(f"  Could not get HF PDB for best variant: {e}. Reusing current baseline.")

            # Elitist selection: champion = highest fitness ever in this lineage; only mutate from champion
            if best_resolved and (champion is None or best_fitness > champion[3]):
                champion = (best_name, best_fasta_dest, best_structure_dest, best_fitness)
                log.info(f"New champion: {best_name} (fitness={best_fitness:.1f})")
            elif champion is not None:
                log.info(f"Keeping champion {champion[0]} (fitness={champion[3]:.1f}); current gen best {best_name} had {best_fitness:.1f}")

            if champion is not None:
                baseline = (champion[0], champion[2], champion[1], crrna_lookup_id)
            else:
                # First gen and best could not be resolved; keep current baseline
                baseline = (baseline_id, baseline_pdb_path, baseline_fasta_path, crrna_lookup_id)

            if gym.mutation_weights:
                bias_file = gym.generate_mpnn_bias_matrix(generation_counter)

        # After MAX_GENERATIONS for this lineage, switch to next
        if not lineage_queue:
            break
        log.info(f"Lineage complete after {MAX_GENERATIONS} generations. Switching to next baseline...")
        baseline, lineage_queue = _get_next_baseline_from_queue(lineage_queue)
        if baseline is None:
            break

    log.info("Evolution loop complete.")

if __name__ == "__main__":
    main_evolution_loop()