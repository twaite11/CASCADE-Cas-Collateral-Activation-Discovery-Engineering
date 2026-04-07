import csv
import os
import json
import multiprocessing
import random
import re
import shutil
import subprocess
import gc
import time
import logging
import importlib.util
import numpy as np

NUM_WORKERS = int(os.environ.get("CASCADE_WORKERS", "3"))


def _flush_gpu_memory():
    """Best-effort GPU VRAM reclaim between heavy inference steps."""
    gc.collect()
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            torch.cuda.synchronize()
    except ImportError:
        pass

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
    EVAL_ENGINE,
    assemble_crrna,
    get_spacer_for_subtype,
)
from utils.pdb_kinematics import calculate_hepn_shift, extract_protenix_scores, find_structure_files


def _find_mini_summary(eval_dir, variant_name, state_suffix):
    """Locate the summary JSON from a mini inference run for a variant's ON or OFF state."""
    import glob as _glob
    pred_dir = os.path.join(eval_dir, f"{variant_name}_{state_suffix}")
    if not os.path.isdir(pred_dir):
        return None, None
    summary_files = _glob.glob(os.path.join(pred_dir, "**", "*_summary*.json"), recursive=True)
    if not summary_files:
        summary_files = _glob.glob(os.path.join(pred_dir, "**", "*_confidence*.json"), recursive=True)
    struct_files = _glob.glob(os.path.join(pred_dir, "**", "*.cif"), recursive=True)
    if not struct_files:
        struct_files = _glob.glob(os.path.join(pred_dir, "**", "*.pdb"), recursive=True)
    if struct_files and summary_files:
        return struct_files[0], summary_files[0]
    return None, None


# --- Configuration ---
METADATA_FILE = "../metadata/variant_domain_metadata.json"
BASE_JSON_DIR = "../jsons"
PHASE1_PDB_DIR = "../outputs/phase1_screening"  # Where Script 2 saved the initial PDBs
GENERATION_DIR = "../outputs/generation_queue"
FAST_EVAL_DIR = "../outputs/fast_eval"
HIGH_FIDELITY_DIR = "../outputs/high_fidelity_scoring"
FINAL_HITS_DIR = "../outputs/optimized_switches"
SWITCH_REPORT_CSV = "../outputs/optimized_switches/switch_report.csv"
GYM_DIR = "../outputs/rl_gym_data"
RL_TRAINING_DATASET = os.path.join(GYM_DIR, "rl_training_dataset.jsonl")
VALIDATED_IDS_FILE = "../outputs/validated_baseline_ids.txt"  # From validate_crispr_repeats.py; restricts lineage to validated repeats
# Optional: set to path to databases/ for Protenix inputprep (improves MSA quality)
SEQRES_DB_PATH = os.path.join(os.path.dirname(__file__), "..", "databases")


_RUST_SEQUTILS_BIN = shutil.which("cascade_sequtils")
if not _RUST_SEQUTILS_BIN:
    _candidate = os.path.join(os.path.dirname(__file__), "..", "rust", "target", "release", "cascade_sequtils")
    if os.name == "nt":
        _candidate += ".exe"
    if os.path.isfile(_candidate):
        _RUST_SEQUTILS_BIN = _candidate


def _run_rust_sequtils(args, timeout=30):
    """Run cascade_sequtils with args. Returns stdout string on success, None on failure."""
    if not _RUST_SEQUTILS_BIN:
        return None
    try:
        r = subprocess.run(
            [_RUST_SEQUTILS_BIN] + args,
            capture_output=True, text=True, timeout=timeout,
        )
        if r.returncode == 0:
            return r.stdout.strip()
    except Exception:
        pass
    return None


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
    # Empty file or validation wrote no IDs (e.g. ViennaRNA missing): do not restrict to zero baselines.
    if not ids:
        log.warning(
            "%s is empty or has no IDs — using all baselines from metadata. "
            "Install ViennaRNA and re-run validation if you need structure-filtered IDs.",
            VALIDATED_IDS_FILE,
        )
        return None
    return ids

# --- Biophysical Thresholds ---
# Distance between catalytic His NE2 atoms (side-chain measurement).
# NE2-NE2 is ~5-7 A shorter than the previous CA-CA metric.
# OFF state: HEPN domains far apart (inactive). ON state: snap together (active).
MIN_OFF_DISTANCE = 18.0  # Ångströms (NE2-NE2; was 25.0 for CA-CA)
MAX_ON_DISTANCE = 12.0   # Ångströms — relaxed from 7.0 to catch near-hits for base-model eval
MIN_IPTM_SCORE = 0.85
MIN_AF2_IG_SCORE = 0.80
# --- Evolution Loop Config ---
MAX_GENERATIONS = 12
VARIANTS_PER_GEN = 5
STAGNATION_LIMIT = 4  # Abandon lineage after this many consecutive gens with no improvement
POPULATION_SIZE = 3   # Top-K population to carry forward across generations
TOURNAMENT_SIZE = 2   # Subset drawn for tournament selection of next parent
MISMATCH_COUNTS = (1, 2, 3)  # Test 1-, 2-, 3-mismatch off-targets; activity at higher count penalized harder
SPECIFICITY_PENALTY_BASE = 0.3  # Base penalty; scaled by mismatch count (3mm > 2mm > 1mm)
# --- Memory / OOM mitigation (seconds; set to 0 to disable) ---
SLEEP_AFTER_PXDESIGN = 2.0       # Allow CUDA driver to reclaim GPU memory after diffusion
SLEEP_AFTER_PROTENIX_BASE = 2.0  # Allow reclaim after heavy base-model ternary prediction
SLEEP_AFTER_PROTENIX_MINI = 0.0  # Optional: use 0.5 if OOM on back-to-back mini runs
FALLBACK_FITNESS_PENALTY = 5.0   # Penalize fallback variants that failed stitching

def compute_fitness(off_dist, on_dist, iptm_score, af2_ig_score, is_full_ternary=False, offtarget_by_mismatch=None):
    """
    Composite fitness for ranking variants.
    Primary metric: non-active when unbound (high OFF dist), highly active when bound (low ON dist, high ipTM).
    AF2-IG: interface confidence for RNA-protein complexes; higher = more reliable ternary prediction.
    Specificity: penalize activity on 1-, 2-, 3-mismatch off-targets; greater mismatch activation = harder penalty.
    """
    multiplier = 2.0 if is_full_ternary else 1.0
    shift = off_dist - on_dist
    fitness = (
        (shift - (MIN_OFF_DISTANCE - MAX_ON_DISTANCE))
        + ((iptm_score - 0.7) * 50)
        + ((af2_ig_score - 0.5) * 20)
    )
    fitness *= multiplier

    # Progressive specificity penalty: activity at 3mm > 2mm > 1mm penalized harder
    if offtarget_by_mismatch:
        for n_mismatch, min_dist in offtarget_by_mismatch.items():
            if min_dist is not None and min_dist < MIN_OFF_DISTANCE:
                penalty = SPECIFICITY_PENALTY_BASE * n_mismatch * (MIN_OFF_DISTANCE - min_dist)
                fitness -= penalty
    return fitness


def _update_population(population, new_results, max_size=POPULATION_SIZE):
    """Merge new generation results into population, keep top-K by fitness.

    Each entry is a dict with keys: name, fasta, fitness, off, on, iptm,
    af2_ig, hf_pdb, crrna_lid, offtarget_by_mismatch.
    """
    combined = list(population) + list(new_results)
    combined.sort(key=lambda x: x["fitness"], reverse=True)
    return combined[:max_size]


def _tournament_select(population, k=TOURNAMENT_SIZE):
    """Pick a parent from population via tournament selection (fitness-proportionate pressure)."""
    if len(population) <= 1:
        return population[0] if population else None
    contestants = random.sample(population, min(k, len(population)))
    return max(contestants, key=lambda x: x["fitness"])


class EvolutionGym:
    """Active Learning Environment for Directed Evolution.

    Uses *relative* fitness: within each generation, the best variant gets
    a positive weight and the worst gets a negative weight.  This ensures
    the bias matrix contains actionable signal even when all absolute
    fitness values are deeply negative (e.g. early exploration).
    """
    def __init__(self, worker_id=0):
        self.worker_id = worker_id
        self.gym_dir = os.path.join(GYM_DIR, f"worker_{worker_id}") if NUM_WORKERS > 1 else GYM_DIR
        os.makedirs(self.gym_dir, exist_ok=True)
        self.mutation_weights = {}
        self.generation_history = []
        self.baseline_fitness = None
        self._pending_evals = []

    def set_baseline_fitness(self, fitness):
        """Set the reference fitness from the unmodified enzyme (Gen 0)."""
        self.baseline_fitness = fitness
        log.info(f"  [RL] Baseline fitness reference set to {fitness:.2f}")

    def register_evaluation(self, variant_id, mutations, off_dist, on_dist, iptm_score, af2_ig_score=0.0, is_full_ternary=False, offtarget_by_mismatch=None):
        """Buffer variant evaluation for end-of-generation relative scoring."""
        fitness = compute_fitness(off_dist, on_dist, iptm_score, af2_ig_score, is_full_ternary, offtarget_by_mismatch)

        self.generation_history.append({
            "variant": variant_id, "fitness": fitness, "mutations": mutations
        })
        self._pending_evals.append({"mutations": mutations, "fitness": fitness})

    def flush_generation(self):
        """Normalize pending evaluations to relative fitness and update weights.

        Relative fitness = (variant_fitness - generation_mean) / spread.
        If a baseline reference exists, center on that instead of the mean.
        This guarantees variants better than average get positive weights
        and worse-than-average get negative — regardless of absolute scale.
        """
        if not self._pending_evals:
            return

        fitnesses = [e["fitness"] for e in self._pending_evals]
        if self.baseline_fitness is not None:
            center = self.baseline_fitness
        else:
            center = sum(fitnesses) / len(fitnesses)

        spread = max(fitnesses) - min(fitnesses) if len(fitnesses) > 1 else 1.0
        spread = max(spread, 1.0)

        for ev in self._pending_evals:
            relative = (ev["fitness"] - center) / spread
            for mut in ev["mutations"]:
                if mut not in self.mutation_weights:
                    self.mutation_weights[mut] = 0.0
                self.mutation_weights[mut] = (self.mutation_weights[mut] * 0.5) + (relative * 0.5)

        n_pos = sum(1 for w in self.mutation_weights.values() if w > 0)
        n_neg = sum(1 for w in self.mutation_weights.values() if w < 0)
        log.info(f"  [RL] Flushed {len(self._pending_evals)} evals → "
                 f"{n_pos} positive / {n_neg} negative mutation weights")
        self._pending_evals = []

    def generate_mpnn_bias_matrix(self, generation_num):
        """Converts weights into a physical bias matrix for PXDesign/ProteinMPNN."""
        bias_matrix = {}
        for mut, weight in self.mutation_weights.items():
            parts = mut.split('_')
            if len(parts) == 2:
                pos = parts[0]
                aa = parts[1]
                if aa != "del" and len(aa) == 1 and aa.isalpha():
                    if pos not in bias_matrix:
                        bias_matrix[pos] = {}
                    bias_matrix[pos][aa] = float(np.clip(weight, -5.0, 5.0))
            
        bias_file = os.path.join(self.gym_dir, f"mpnn_bias_gen_{generation_num}.json")
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
    structure_path, offtarget_by_mismatch, is_elite, rl_dataset_path=None,
):
    """
    Append a single variant evaluation to rl_training_dataset.jsonl.
    Format is designed for DRAKES/ProteinMPNN post-training: (structure, sequence, reward).
    """
    dest = rl_dataset_path or RL_TRAINING_DATASET
    os.makedirs(os.path.dirname(dest), exist_ok=True)
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
    with open(dest, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")


def build_metadata_override_for_evolved(baseline_id, baseline_fasta_path, crrna_lookup_id, domain_metadata):
    """Build metadata override dict for evolved variants not in variant_domain_metadata.json.
    Uses _select_hepn_pair for consistent HEPN domain pairing (not first/last regex hit)."""
    with open(baseline_fasta_path, 'r') as f:
        seq = "".join([l.strip() for l in f if not l.startswith(">")])
    motif = re.compile(r'R.{3,6}H')
    pair = _select_hepn_pair(seq, motif)
    if pair is None:
        return None
    hepn1_center = pair[0].start()
    hepn2_center = pair[1].start()
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
            "subtype": parent_data.get("subtype", "unknown"),
        }
    }


def save_crrna_for_elite(variant_name, crrna_lookup_id, domain_metadata):
    """Saves the crRNA sequence (repeat + spacer) for an elite switch to FINAL_HITS_DIR."""
    parent = domain_metadata.get(crrna_lookup_id)
    if not parent:
        return
    subtype = parent.get("subtype", "unknown")
    spacer = get_spacer_for_subtype(subtype)
    crrna_seq = assemble_crrna(parent["crRNA_repeat_used"], spacer, subtype)
    os.makedirs(FINAL_HITS_DIR, exist_ok=True)
    crrna_path = os.path.join(FINAL_HITS_DIR, f"{variant_name}_crRNA.fasta")
    with open(crrna_path, 'w') as f:
        f.write(f">{variant_name}_crRNA\n")
        f.write(f"{crrna_seq}\n")


_REPORT_FIELDS = [
    "variant", "lineage", "generation", "off_dist_A", "on_dist_A", "delta_A",
    "iptm", "af2_ig", "fitness", "n_mutations", "filter_passed", "is_elite",
    "score_source", "n_hepn_motifs", "hepn_spacing", "subtype",
]
_report_lock = None


def _init_switch_report():
    """Create the CSV header if the report doesn't exist yet."""
    report_path = os.path.join(os.path.dirname(__file__), SWITCH_REPORT_CSV)
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    if not os.path.exists(report_path):
        with open(report_path, "w", newline="", encoding="utf-8") as f:
            csv.writer(f).writerow(_REPORT_FIELDS)


def append_switch_report(variant_name, lineage_id, generation, off_dist, on_dist,
                         iptm, af2_ig, fitness, n_mutations, filter_passed,
                         is_elite, score_source, protein_seq=None, subtype="unknown"):
    """Append one row to the switch report CSV (thread-safe via lock)."""
    n_hepn = ""
    hepn_spacing = ""
    if protein_seq:
        motif = re.compile(r'R.{3,6}H')
        matches = list(motif.finditer(protein_seq))
        n_hepn = len(matches)
        pair = _select_hepn_pair(protein_seq, motif)
        if pair:
            hepn_spacing = pair[1].start() - pair[0].start()

    row = [
        variant_name, lineage_id, generation,
        f"{off_dist:.1f}", f"{on_dist:.1f}", f"{off_dist - on_dist:.1f}",
        f"{iptm:.4f}", f"{af2_ig:.4f}", f"{fitness:.2f}",
        n_mutations, filter_passed, is_elite, score_source,
        n_hepn, hepn_spacing, subtype,
    ]

    report_path = os.path.join(os.path.dirname(__file__), SWITCH_REPORT_CSV)
    try:
        with open(report_path, "a", newline="", encoding="utf-8") as f:
            csv.writer(f).writerow(row)
    except OSError:
        pass


def _select_hepn_pair(sequence, motif, min_sep=150, max_sep=600, ideal_sep=300):
    """
    Select the best pair of R...H HEPN catalytic motifs from a protein sequence.

    Uses positional and separation constraints rather than naive first/last:
      - HEPN1 expected in the first 65% of the protein
      - HEPN2 expected in the last 65%
      - Typical separation: 150-600 residues (~300 ideal)

    Returns (match1, match2) regex match objects, or None.
    """
    matches = list(motif.finditer(sequence))
    if len(matches) < 2:
        return None

    seq_len = len(sequence)
    early = [m for m in matches if m.start() < seq_len * 0.65]
    late = [m for m in matches if m.start() > seq_len * 0.35]

    best_pair = None
    best_score = float('inf')
    for m1 in early:
        for m2 in late:
            sep = m2.start() - m1.start()
            if min_sep <= sep <= max_sep:
                score = abs(sep - ideal_sep)
                if score < best_score:
                    best_score = score
                    best_pair = (m1, m2)

    if best_pair is not None:
        return best_pair

    for i, m1 in enumerate(matches):
        for m2 in matches[i + 1:]:
            sep = m2.start() - m1.start()
            if min_sep <= sep <= max_sep:
                score = abs(sep - ideal_sep)
                if score < best_score:
                    best_score = score
                    best_pair = (m1, m2)

    return best_pair


def get_catalytic_histidine_indices(fasta_path):
    """Parses a FASTA to find the exact 1-based indices of the two catalytic Histidines.
    Uses only the first sequence if the FASTA contains multiple entries.
    Always uses _select_hepn_pair for consistent pairing with 01_parse_and_annotate."""
    with open(fasta_path, 'r') as f:
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_lines:
                    break
                continue
            seq_lines.append(line)
        seq = "".join(seq_lines)

    motif = re.compile(r'R.{3,6}H')
    pair = _select_hepn_pair(seq, motif)
    if pair is None:
        return None, None

    return pair[0].end(), pair[1].end()

def extract_mutations(baseline_id, variant_fasta, baseline_fasta_path=None):
    """Compares the new variant against the original sequence to map the mutations.
    Handles length mismatches (indels) by padding the shorter sequence.
    baseline_fasta_path: if provided, load baseline sequence from FASTA (for evolved baselines)."""
    baseline_file = baseline_fasta_path if (baseline_fasta_path and os.path.exists(baseline_fasta_path)) else os.path.join(BASE_JSON_DIR, f"{baseline_id}.json")
    rust_out = _run_rust_sequtils(["extract-mutations", "--baseline", str(baseline_file), "--variant", str(variant_fasta)])
    if rust_out is not None:
        try:
            return json.loads(rust_out)
        except (json.JSONDecodeError, TypeError):
            pass

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

def _find_existing_structure(final_dir, name):
    """Search for a ternary complex structure (CIF preferred, then PDB)."""
    base = os.path.join(final_dir, f"{name}_ternary_complex")
    for ext in (".cif", ".pdb"):
        p = base + ext
        if os.path.exists(p):
            return p
    return None


def _resolve_best_baseline(best_name, best_fasta, best_hf_pdb, crrna_lookup_id, domain_metadata):
    """
    Resolve the best variant into a usable baseline tuple:
    (baseline_id, pdb_path, fasta_path, crrna_lookup_id).
    Copies files to FINAL_HITS_DIR and returns the new baseline or None on failure.
    """
    if not best_fasta or not os.path.exists(best_fasta):
        log.warning(f"No FASTA available for {best_name}; cannot resolve as baseline")
        return None
    os.makedirs(FINAL_HITS_DIR, exist_ok=True)
    best_fasta_dest = os.path.join(FINAL_HITS_DIR, f"{best_name}_optimal.fasta")
    shutil.copy(best_fasta, best_fasta_dest)
    save_crrna_for_elite(best_name, crrna_lookup_id, domain_metadata)

    if best_hf_pdb:
        ext = os.path.splitext(best_hf_pdb)[1] or ".pdb"
        best_structure_dest = os.path.join(FINAL_HITS_DIR, f"{best_name}_ternary_complex{ext}")
        shutil.copy(best_hf_pdb, best_structure_dest)
        return (best_name, best_structure_dest, best_fasta_dest, crrna_lookup_id)

    existing = _find_existing_structure(FINAL_HITS_DIR, best_name)
    if existing:
        return (best_name, existing, best_fasta_dest, crrna_lookup_id)

    # Try to generate a high-fidelity structure for this variant
    on_json = os.path.join(FAST_EVAL_DIR, f"{best_name}_ON.json")
    if os.path.exists(on_json):
        try:
            hf_structure, _ = run_protenix_inference(
                on_json, HIGH_FIDELITY_DIR, model_tier="base", seqres_db_path=SEQRES_DB_PATH
            )
            ext = os.path.splitext(hf_structure)[1] or ".pdb"
            hf_dest = os.path.join(FINAL_HITS_DIR, f"{best_name}_ternary_complex{ext}")
            shutil.copy(hf_structure, hf_dest)
            return (best_name, hf_dest, best_fasta_dest, crrna_lookup_id)
        except Exception as e:
            log.warning(f"Could not get HF PDB for best variant {best_name}: {e}")

    return None


def _evaluate_baseline_reference(baseline_id, baseline_fasta_path, crrna_lookup_id,
                                  domain_metadata, gpu_lock, fast_eval_dir):
    """
    Evaluate the unmodified baseline enzyme (Gen 0) to establish a reference
    fitness for the RL system.
    """
    log.info("=" * 60)
    log.info(f"Gen 0 — Evaluating unmodified baseline: {baseline_id}")
    log.info("=" * 60)

    seq = _read_baseline_sequence(baseline_id, baseline_fasta_path)
    if not seq:
        log.warning(f"Could not read baseline sequence for {baseline_id}")
        return None

    motif_check = re.compile(r'R.{3,6}H')
    n_motifs = len(list(motif_check.finditer(seq)))
    if n_motifs > 8:
        log.warning(f"Baseline {baseline_id} has {n_motifs} R...H motifs — skipping (likely not a clean Cas13)")
        return None
    if n_motifs > 4:
        log.info(f"Baseline {baseline_id}: {n_motifs} R...H motifs (expected 2 for Cas13 — selecting best pair)")

    import tempfile
    baseline_fasta_dir = os.path.join(FAST_EVAL_DIR, "baseline_fastas")
    os.makedirs(baseline_fasta_dir, exist_ok=True)
    persistent_fasta = os.path.join(baseline_fasta_dir, f"{baseline_id}.fasta")
    with open(persistent_fasta, "w") as f:
        f.write(f">{baseline_id}\n{seq}\n")
    tmp_fasta = persistent_fasta

    try:
        h1_idx, h2_idx = get_catalytic_histidine_indices(tmp_fasta)
        if h1_idx is None or h2_idx is None:
            log.warning(f"Could not find HEPN motifs in baseline {baseline_id}")
            return None

        off_json, on_json = generate_evaluation_jsons(
            tmp_fasta, baseline_id, METADATA_FILE, fast_eval_dir, crrna_lookup_id=crrna_lookup_id
        )

        with gpu_lock:
            off_pdb, off_summary = run_protenix_inference(
                off_json, fast_eval_dir, model_tier="mini", seqres_db_path=SEQRES_DB_PATH
            )
        with gpu_lock:
            on_pdb, on_summary = run_protenix_inference(
                on_json, fast_eval_dir, model_tier="mini", seqres_db_path=SEQRES_DB_PATH
            )

        off_dist = calculate_hepn_shift(off_pdb, h1_idx, h2_idx)
        on_dist = calculate_hepn_shift(on_pdb, h1_idx, h2_idx)

        if on_summary:
            scores = extract_protenix_scores(on_summary)
            iptm = scores["iptm"] if scores["iptm"] > 0.0 else 0.4
            af2_ig = scores["af2_ig"]
        else:
            iptm, af2_ig = 0.4, 0.0

        fitness = compute_fitness(off_dist, on_dist, iptm, af2_ig, False, None)
        log.info(
            f"[Gen 0] Baseline {baseline_id}: OFF={off_dist:.1f}A ON={on_dist:.1f}A "
            f"delta={off_dist - on_dist:.1f}A iptm={iptm:.3f} fitness={fitness:.2f}"
        )
        return {
            "name": baseline_id, "fasta": persistent_fasta, "fitness": fitness,
            "off": off_dist, "on": on_dist, "iptm": iptm, "af2_ig": af2_ig,
            "hf_pdb": None, "crrna_lid": crrna_lookup_id, "offtarget": None,
        }

    except Exception as e:
        log.warning(f"Baseline evaluation failed for {baseline_id}: {e}")
        return None


def _run_single_lineage(worker_id, baselines, gpu_lock):
    """
    Run one or more lineages end-to-end. Fully self-contained with its own
    EvolutionGym and output directories. The gpu_lock serializes Protenix
    inference across parallel workers.

    baselines: list of (baseline_id, pdb_path, fasta_path, crrna_lookup_id) tuples.
    """
    tag = f"[W{worker_id}]"
    log.info(f"{tag} Starting worker with {len(baselines)} lineage(s)")

    with open(METADATA_FILE, 'r') as f:
        domain_metadata = json.load(f)

    gym = EvolutionGym(worker_id=worker_id)

    suffix = f"worker_{worker_id}" if NUM_WORKERS > 1 else ""
    fast_eval_dir = os.path.join(FAST_EVAL_DIR, suffix) if suffix else FAST_EVAL_DIR
    hf_dir = os.path.join(HIGH_FIDELITY_DIR, suffix) if suffix else HIGH_FIDELITY_DIR
    gen_dir = os.path.join(GENERATION_DIR, suffix) if suffix else GENERATION_DIR
    rl_dataset = os.path.join(gym.gym_dir, "rl_training_dataset.jsonl")
    os.makedirs(fast_eval_dir, exist_ok=True)
    os.makedirs(hf_dir, exist_ok=True)
    os.makedirs(gen_dir, exist_ok=True)

    mismatch_seqs = generate_mismatch_sequences(TARGET_REGION, mismatch_counts=MISMATCH_COUNTS, num_per_count=3, seed=42)

    global_best = None
    bias_file = None

    for baseline in baselines:
        baseline_id, baseline_pdb_path, baseline_fasta_path, crrna_lookup_id = baseline

        # Gen 0: evaluate the unmodified baseline
        baseline_result = _evaluate_baseline_reference(
            baseline_id, baseline_fasta_path, crrna_lookup_id, domain_metadata, gpu_lock, fast_eval_dir,
        )
        if baseline_result is not None:
            gym.set_baseline_fitness(baseline_result["fitness"])

        current_baseline = baseline
        # Seed population with the baseline so variants must beat it to survive
        population = [baseline_result] if baseline_result is not None else []
        stagnation_counter = 0
        lineage_elite_found = False

        for generation_counter in range(1, MAX_GENERATIONS + 1):
            # Tournament-select a parent from the population (falls back to current_baseline if empty)
            if population:
                parent = _tournament_select(population, k=TOURNAMENT_SIZE)
                bid = parent["name"]
                bpdb = parent.get("hf_pdb") or current_baseline[1]
                bfasta = parent["fasta"]
                crrna_lid = parent["crrna_lid"]
            else:
                bid, bpdb, bfasta, crrna_lid = current_baseline

            log.info("=" * 60)
            log.info(f"{tag} Generation {generation_counter} | Parent: {bid}")
            pop_summary = ", ".join(f"{p['name']}={p['fitness']:.1f}" for p in population[:3])
            log.info(f"{tag}   Population[{len(population)}]: [{pop_summary}]")
            if global_best:
                log.info(f"{tag}   Global best: {global_best[0]} (fitness={global_best[2]:.1f})")
            log.info("=" * 60)

            metadata_override = None
            if bfasta:
                metadata_override = build_metadata_override_for_evolved(bid, bfasta, crrna_lid, domain_metadata)
                if not metadata_override:
                    log.warning(f"{tag} Could not build metadata for evolved baseline {bid}. Skipping lineage.")
                    break

            try:
                new_variants_fastas = pxdesign_wrapper.run_pxdesign_generation(
                    baseline_structure=bpdb,
                    variant_id=bid,
                    metadata_path=METADATA_FILE,
                    bias_json_path=bias_file,
                    output_dir=os.path.join(gen_dir, f"gen_{generation_counter}"),
                    variant_count=VARIANTS_PER_GEN,
                    metadata_override=metadata_override,
                    baseline_fasta_path=bfasta,
                    base_json_dir=BASE_JSON_DIR,
                    generation_num=generation_counter,
                    lineage_seed=crrna_lid,
                )
            except Exception as e:
                log.error(f"{tag} PXDesign failed: {e}. Skipping generation...")
                continue

            if SLEEP_AFTER_PXDESIGN > 0:
                time.sleep(SLEEP_AFTER_PXDESIGN)
            _flush_gpu_memory()

            if not new_variants_fastas:
                log.warning(f"{tag} No variants generated. Continuing...")
                continue

            results = []

            for variant_fasta in new_variants_fastas:
                mutations_made = extract_mutations(bid, variant_fasta, bfasta)
                variant_name = os.path.basename(variant_fasta).replace(".fasta", "")

                h1_idx, h2_idx = get_catalytic_histidine_indices(variant_fasta)
                if h1_idx is None or h2_idx is None:
                    continue

                log.info(f"{tag} Evaluating {variant_name} (mini OFF/ON)...")
                off_json, on_json = generate_evaluation_jsons(
                    variant_fasta, bid, METADATA_FILE, fast_eval_dir, crrna_lookup_id=crrna_lid
                )

                try:
                    with gpu_lock:
                        off_pdb, off_summary = run_protenix_inference(
                            off_json, fast_eval_dir, model_tier="mini", seqres_db_path=SEQRES_DB_PATH
                        )
                    if SLEEP_AFTER_PROTENIX_MINI > 0:
                        time.sleep(SLEEP_AFTER_PROTENIX_MINI)
                    with gpu_lock:
                        on_pdb, on_summary = run_protenix_inference(
                            on_json, fast_eval_dir, model_tier="mini", seqres_db_path=SEQRES_DB_PATH
                        )
                    if SLEEP_AFTER_PROTENIX_MINI > 0:
                        time.sleep(SLEEP_AFTER_PROTENIX_MINI)
                except Exception as e:
                    log.warning(f"{tag} {EVAL_ENGINE} failed for {variant_name}: {e}")
                    fitness = compute_fitness(0, 999, 0.4, 0, False, None)
                    gym.register_evaluation(variant_name, mutations_made, 0, 999, 0.4, af2_ig_score=0.0, is_full_ternary=False)
                    save_rl_training_record(
                        variant_name, variant_fasta, bid, bfasta, crrna_lid,
                        generation_counter, mutations_made, fitness, 0, 999, 0.4, 0.0,
                        None, None, False, rl_dataset_path=rl_dataset,
                    )
                    results.append({"name": variant_name, "fasta": variant_fasta, "fitness": fitness,
                                    "off": 0, "on": 999, "iptm": 0.4, "af2_ig": 0.0,
                                    "hf_pdb": None, "crrna_lid": crrna_lid, "offtarget": None})
                    continue

                off_dist = calculate_hepn_shift(off_pdb, h1_idx, h2_idx)
                on_dist = calculate_hepn_shift(on_pdb, h1_idx, h2_idx)

                has_potential = (off_dist >= MIN_OFF_DISTANCE) and (on_dist <= MAX_ON_DISTANCE)
                log.info(f"{tag} [HEPN mini] {variant_name} OFF={off_dist:.1f}A ON={on_dist:.1f}A delta={off_dist - on_dist:.1f}A")
                log.info(f"{tag} [FilterGate] {variant_name} pass={has_potential} (OFF>={MIN_OFF_DISTANCE}A and ON<={MAX_ON_DISTANCE}A)")

                offtarget_by_mismatch = {}
                hf_pdb_path = None
                true_on_dist = on_dist

                if on_summary:
                    mini_scores = extract_protenix_scores(on_summary)
                    iptm = mini_scores["iptm"] if mini_scores["iptm"] > 0.0 else 0.4
                    af2_ig = mini_scores["af2_ig"]
                else:
                    iptm, af2_ig = 0.4, 0.0

                if has_potential:
                    log.info(f"{tag} Filter passed. Running base ternary...")
                    with gpu_lock:
                        hf_pdb, hf_summary = run_protenix_inference(
                            on_json, hf_dir, model_tier="base", seqres_db_path=SEQRES_DB_PATH
                        )
                    if SLEEP_AFTER_PROTENIX_BASE > 0:
                        time.sleep(SLEEP_AFTER_PROTENIX_BASE)
                    _flush_gpu_memory()
                    true_on_dist = calculate_hepn_shift(hf_pdb, h1_idx, h2_idx)
                    scores = extract_protenix_scores(hf_summary)
                    iptm, af2_ig = scores["iptm"], scores["af2_ig"]
                    hf_pdb_path = hf_pdb

                    for i, (ot_rna, n_mismatch) in enumerate(mismatch_seqs):
                        try:
                            ot_json = generate_offtarget_json(
                                variant_fasta, crrna_lid, METADATA_FILE, ot_rna, fast_eval_dir,
                                suffix=f"{n_mismatch}mm_{i}"
                            )
                            with gpu_lock:
                                ot_pdb, _ = run_protenix_inference(
                                    ot_json, fast_eval_dir, model_tier="mini", seqres_db_path=SEQRES_DB_PATH
                                )
                            if SLEEP_AFTER_PROTENIX_MINI > 0:
                                time.sleep(SLEEP_AFTER_PROTENIX_MINI)
                            ot_dist = calculate_hepn_shift(ot_pdb, h1_idx, h2_idx)
                            if n_mismatch not in offtarget_by_mismatch:
                                offtarget_by_mismatch[n_mismatch] = ot_dist
                            else:
                                offtarget_by_mismatch[n_mismatch] = min(offtarget_by_mismatch[n_mismatch], ot_dist)
                        except Exception:
                            offtarget_by_mismatch[n_mismatch] = MIN_OFF_DISTANCE
                    if offtarget_by_mismatch:
                        mm_str = " | ".join(f"{k}mm:{v:.1f}A" for k, v in sorted(offtarget_by_mismatch.items()))
                        log.info(f"{tag} [Specificity] {variant_name} | {mm_str}")

                score_source = "base" if has_potential else ("mini" if on_summary else "default")
                fitness = compute_fitness(off_dist, true_on_dist, iptm, af2_ig, has_potential, offtarget_by_mismatch or None)
                log.info(
                    f"{tag} [HEPN scored] {variant_name} OFF={off_dist:.1f}A ON={true_on_dist:.1f}A "
                    f"delta={off_dist - true_on_dist:.1f}A iptm={iptm:.3f} af2_ig={af2_ig:.3f} "
                    f"fitness={fitness:.2f} (scores from {score_source})"
                )

                if "fallback" in variant_name:
                    fitness -= FALLBACK_FITNESS_PENALTY
                gym.register_evaluation(
                    variant_name, mutations_made, off_dist, true_on_dist, iptm,
                    af2_ig_score=af2_ig, is_full_ternary=has_potential,
                    offtarget_by_mismatch=offtarget_by_mismatch or None
                )
                struct_path = hf_pdb_path if hf_pdb_path else on_pdb
                is_elite = (iptm >= MIN_IPTM_SCORE and af2_ig >= MIN_AF2_IG_SCORE and true_on_dist <= MAX_ON_DISTANCE)
                save_rl_training_record(
                    variant_name, variant_fasta, bid, bfasta, crrna_lid,
                    generation_counter, mutations_made, fitness, off_dist, true_on_dist, iptm, af2_ig,
                    struct_path, offtarget_by_mismatch or None, is_elite, rl_dataset_path=rl_dataset,
                )

                variant_protein = None
                try:
                    with open(variant_fasta, 'r') as _vf:
                        variant_protein = "".join(l.strip() for l in _vf if not l.startswith(">"))
                except OSError:
                    pass
                _sub = domain_metadata.get(crrna_lid, {}).get("subtype", "unknown")
                append_switch_report(
                    variant_name, bid, generation_counter, off_dist, true_on_dist,
                    iptm, af2_ig, fitness, len(mutations_made) if mutations_made else 0,
                    has_potential, is_elite, score_source,
                    protein_seq=variant_protein, subtype=_sub,
                )

                results.append({"name": variant_name, "fasta": variant_fasta, "fitness": fitness,
                                "off": off_dist, "on": true_on_dist, "iptm": iptm, "af2_ig": af2_ig,
                                "hf_pdb": hf_pdb_path, "crrna_lid": crrna_lid, "offtarget": offtarget_by_mismatch})
                _flush_gpu_memory()

            if not results:
                log.warning(f"{tag} No valid results this generation. Continuing...")
                continue

            gen_best = max(results, key=lambda r: r["fitness"])
            log.info(f"{tag} Generation best: {gen_best['name']} (fitness={gen_best['fitness']:.1f})")

            # --- Population update: merge new results into top-K pool ---
            prev_best_fitness = population[0]["fitness"] if population else None
            population = _update_population(population, results, max_size=POPULATION_SIZE)
            log.info(f"{tag} Population after merge: "
                     + ", ".join(f"{p['name']}={p['fitness']:.1f}" for p in population))

            if gen_best["iptm"] >= MIN_IPTM_SCORE and gen_best["af2_ig"] >= MIN_AF2_IG_SCORE and gen_best["on"] <= MAX_ON_DISTANCE:
                log.info(f"{tag} ELITE TERNARY SWITCH FOUND!")
                os.makedirs(FINAL_HITS_DIR, exist_ok=True)
                shutil.copy(gen_best["fasta"], os.path.join(FINAL_HITS_DIR, f"{gen_best['name']}_optimal.fasta"))
                if gen_best["hf_pdb"]:
                    ext = os.path.splitext(gen_best["hf_pdb"])[1] or ".pdb"
                    shutil.copy(gen_best["hf_pdb"], os.path.join(FINAL_HITS_DIR, f"{gen_best['name']}_ternary_complex{ext}"))
                save_crrna_for_elite(gen_best["name"], crrna_lid, domain_metadata)
                lineage_elite_found = True

            # Update global best from the population leader
            pop_leader = population[0]
            if global_best is None or pop_leader["fitness"] > global_best[2]:
                resolved = _resolve_best_baseline(
                    pop_leader["name"], pop_leader["fasta"], pop_leader["hf_pdb"], crrna_lid, domain_metadata
                )
                if resolved:
                    global_best = (
                        pop_leader["name"], resolved[2], pop_leader["fitness"],
                        pop_leader["off"], pop_leader["on"], pop_leader["iptm"], pop_leader["af2_ig"],
                        resolved[1], crrna_lid,
                    )
                    current_baseline = resolved
                    log.info(f"{tag} NEW GLOBAL BEST: {pop_leader['name']} (fitness={pop_leader['fitness']:.1f})")

            # Stagnation: check if the population leader improved
            new_best_fitness = population[0]["fitness"]
            if prev_best_fitness is not None and new_best_fitness <= prev_best_fitness:
                stagnation_counter += 1
            else:
                stagnation_counter = 0

            gym.flush_generation()
            if gym.mutation_weights:
                bias_file = gym.generate_mpnn_bias_matrix(generation_counter)

            if lineage_elite_found:
                log.info(f"{tag} Elite found — moving to next lineage")
                break

            if stagnation_counter >= STAGNATION_LIMIT:
                log.info(f"{tag} No improvement for {STAGNATION_LIMIT} consecutive generations — abandoning lineage")
                break

        log.info(f"{tag} Lineage {baseline_id} complete after {generation_counter} generation(s).")

    log.info(f"{tag} Worker finished all {len(baselines)} lineage(s).")
    return global_best


def _aggregate_worker_results():
    """Merge per-worker RL training datasets into one file."""
    merged_path = RL_TRAINING_DATASET
    os.makedirs(os.path.dirname(merged_path), exist_ok=True)
    import glob as _glob
    worker_files = sorted(_glob.glob(os.path.join(GYM_DIR, "worker_*", "rl_training_dataset.jsonl")))
    if not worker_files:
        return
    with open(merged_path, "a", encoding="utf-8") as out:
        for wf in worker_files:
            with open(wf, encoding="utf-8") as inp:
                for line in inp:
                    out.write(line)
    log.info(f"Merged {len(worker_files)} worker RL datasets into {merged_path}")


class _DummyLock:
    """No-op context manager for single-worker mode (avoids multiprocessing overhead)."""
    def __enter__(self):
        return self
    def __exit__(self, *args):
        pass


def main_evolution_loop():
    log.info("Initializing SwitchBlade Active Learning Evolution Loop...")
    log.info(f"Workers: {NUM_WORKERS} (set CASCADE_WORKERS env to change)")
    os.makedirs(FAST_EVAL_DIR, exist_ok=True)
    os.makedirs(HIGH_FIDELITY_DIR, exist_ok=True)
    os.makedirs(GYM_DIR, exist_ok=True)
    _init_switch_report()
    log.info(f"RL training data will be appended to {RL_TRAINING_DATASET} (see RL_TRAINING_FORMAT.md)")

    with open(METADATA_FILE, 'r') as f:
        domain_metadata = json.load(f)

    baseline_ids = list(domain_metadata.keys())
    validated = _load_validated_baseline_ids()
    if validated is not None:
        baseline_ids = [b for b in baseline_ids if b in validated]
        log.info(f"Restricting to {len(baseline_ids)} baselines with validated CRISPR repeats (from {VALIDATED_IDS_FILE})")
    else:
        log.info(f"Using all {len(baseline_ids)} baselines (not restricting by {VALIDATED_IDS_FILE})")

    raw_queue = [(bid, None, None, bid) for bid in baseline_ids]
    log.info(f"Queue contains {len(raw_queue)} lineages ({MAX_GENERATIONS} generations each)")
    if not raw_queue:
        log.warning("No baselines in metadata. Exiting.")
        return

    resolved_baselines = []
    remaining = list(raw_queue)
    while remaining:
        baseline, remaining = _get_next_baseline_from_queue(remaining)
        if baseline is None:
            break
        resolved_baselines.append(baseline)

    if not resolved_baselines:
        log.warning("No valid Phase 1 structures found for any baseline. Exiting.")
        return

    log.info(f"Resolved {len(resolved_baselines)} baselines with Phase 1 structures")

    n_workers = min(NUM_WORKERS, len(resolved_baselines))

    if n_workers <= 1:
        log.info("Running single-worker mode")
        gpu_lock = _DummyLock()
        _run_single_lineage(0, resolved_baselines, gpu_lock)
    else:
        chunks = [resolved_baselines[i::n_workers] for i in range(n_workers)]
        log.info(f"Distributing {len(resolved_baselines)} lineages across {n_workers} workers: "
                 + ", ".join(f"W{i}={len(c)}" for i, c in enumerate(chunks)))

        gpu_lock = multiprocessing.Manager().Lock()

        with multiprocessing.Pool(n_workers) as pool:
            pool.starmap(_run_single_lineage, [
                (i, chunks[i], gpu_lock) for i in range(n_workers)
            ])

        _aggregate_worker_results()

    log.info("Evolution loop complete.")


if __name__ == "__main__":
    main_evolution_loop()