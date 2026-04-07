import csv
import json
import os
import re
import glob
import shutil
import sqlite3
import subprocess
import sys
import logging

csv.field_size_limit(sys.maxsize)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# --- Configuration ---
DATA_DIR = "../data/mined_hits"
DB_FILE = "../metadata/cas13_variants.db"
JSON_OUT_DIR = "../jsons"
METADATA_OUT_FILE = "../metadata/variant_domain_metadata.json"

# RNA Constants — use the subtype-aware helpers from protenix_eval when available,
# but keep these defaults for standalone use.
DUMMY_SPACER_RNA = "GUCGACUGACGUACGUACGUACGU" # 24nt fallback
DUMMY_TARGET_RNA = "AAAAAA" + "ACGUACGUACGUACGUCAGUCGAC" + "AAAAAA"

try:
    from utils.protenix_eval import (
        assemble_crrna as _assemble_crrna,
        get_spacer_for_subtype as _get_spacer,
        get_target_for_spacer as _get_target,
    )
    _HAS_CRRNA_HELPERS = True
except ImportError:
    _HAS_CRRNA_HELPERS = False
BATCH_SIZE = 1000  # Number of sequences to hold in memory at once

def init_db():
    """Initializes the local SQLite database schema."""
    os.makedirs(os.path.dirname(DB_FILE), exist_ok=True)
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS variants (
            sequence_id TEXT PRIMARY KEY,
            sra_accession TEXT,
            sequence TEXT,
            crrna_repeat TEXT,
            score REAL,
            hepn1_start INTEGER,
            hepn1_end INTEGER,
            hepn2_start INTEGER,
            hepn2_end INTEGER,
            status TEXT,
            reason TEXT
        )
    ''')
    conn.commit()
    return conn

def load_files_to_db(conn):
    """Parses all FASTA and CSV files in the DATA_DIR and streams them into SQLite."""
    if not os.path.exists(DATA_DIR):
        log.warning(f"Directory {DATA_DIR} does not exist. Creating it now.")
        os.makedirs(DATA_DIR, exist_ok=True)
        return 0
        
    fasta_files = glob.glob(os.path.join(DATA_DIR, "*.fasta"))
    csv_files = glob.glob(os.path.join(DATA_DIR, "*.csv"))
    
    cursor = conn.cursor()
    items_processed = 0
    
    # 1. Parse all FASTAs and stream base sequences
    log.info("Streaming FASTA sequences to SQLite...")
    for fasta in fasta_files:
        current_id = ""
        current_seq = []
        with open(fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id:
                        # Upsert logic to prevent duplicate crashes if run multiple times
                        cursor.execute('''
                            INSERT INTO variants (sequence_id, sequence) 
                            VALUES (?, ?) 
                            ON CONFLICT(sequence_id) DO UPDATE SET sequence=excluded.sequence
                        ''', (current_id, "".join(current_seq)))
                        items_processed += 1
                    current_id = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            # Catch the last sequence
            if current_id:
                cursor.execute('''
                    INSERT INTO variants (sequence_id, sequence) 
                    VALUES (?, ?) 
                    ON CONFLICT(sequence_id) DO UPDATE SET sequence=excluded.sequence
                ''', (current_id, "".join(current_seq)))
                items_processed += 1
    conn.commit()
                
    # 2. Parse all CSVs and update metadata
    log.info("Streaming CSV metadata to SQLite...")

    # If fix_crrna_assignments has produced a corrected report, use it
    corrected_dr_path = os.path.join(os.path.dirname(__file__), "..", "outputs", "repeat_validation_report.csv")
    corrected_drs = {}
    if os.path.exists(corrected_dr_path):
        try:
            with open(corrected_dr_path, encoding="utf-8") as cf:
                cr = csv.DictReader(cf)
                for crow in cr:
                    sid = crow.get("sequence_id", "")
                    new_dr = crow.get("new_dr", "")
                    if sid and new_dr:
                        corrected_drs[sid] = new_dr
            if corrected_drs:
                log.info(f"Loaded {len(corrected_drs)} corrected crRNA DRs from {corrected_dr_path}")
        except Exception:
            pass

    for csv_file in csv_files:
        with open(csv_file, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                seq_id = row.get("sequence_id")
                sra_acc = row.get("sra_accession", "")
                score = float(row.get("score", 0.0))

                if seq_id in corrected_drs:
                    repeat_str = corrected_drs[seq_id]
                else:
                    repeat_str = row.get("repeat_domains", "").split("|")[0]

                cursor.execute('''
                    UPDATE variants 
                    SET crrna_repeat = ?, sra_accession = ?, score = ?
                    WHERE sequence_id = ?
                ''', (repeat_str, sra_acc, score, seq_id))
    conn.commit()
                
    return items_processed

def _select_hepn_pair(sequence, motif, min_sep=150, max_sep=600, ideal_sep=300):
    """
    Select the best pair of R...H HEPN catalytic motifs from a protein sequence.

    Instead of naively taking first/last regex hit (which anchors on spurious
    R...H occurrences in NTD or C-terminal regions), this applies positional
    and separation constraints based on known Cas13 domain architecture:
      - HEPN1 is expected in the first 65% of the protein
      - HEPN2 is expected in the last 65% of the protein
      - The two domains are typically 150-600 residues apart (~300 typical)

    Returns (hepn1_center, hepn2_center) as 0-based positions, or None.
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
                    best_pair = (m1.start(), m2.start())

    if best_pair is not None:
        return best_pair

    # Fallback: if no pair satisfies the strict positional constraints,
    # try all pairs with just the separation filter
    for i, m1 in enumerate(matches):
        for m2 in matches[i + 1:]:
            sep = m2.start() - m1.start()
            if min_sep <= sep <= max_sep:
                score = abs(sep - ideal_sep)
                if score < best_score:
                    best_score = score
                    best_pair = (m1.start(), m2.start())

    return best_pair


def identify_hepn_domains(conn):
    """
    Scans sequences in the DB strictly for HEPN domains. 
    Processes sequentially using fetchmany() to prevent memory overloading.
    """
    log.info("Scanning sequences for HEPN domains (may take a minute for large DBs)...")
    
    cursor = conn.cursor()
    update_cursor = conn.cursor()
    
    # Select only sequences that haven't been successfully processed yet
    cursor.execute("SELECT sequence_id, sequence FROM variants WHERE sequence IS NOT NULL")
    
    motif = re.compile(r'R.{3,6}H')
    processed = 0
    
    while True:
        batch = cursor.fetchmany(BATCH_SIZE)
        if not batch:
            break
            
        for seq_id, sequence in batch:
            pair = _select_hepn_pair(sequence, motif)
            
            if pair is None:
                matches = list(motif.finditer(sequence))
                update_cursor.execute('''
                    UPDATE variants SET status = 'failed', reason = ? WHERE sequence_id = ?
                ''', (f"Only {len(matches)} HEPN motifs found (or no valid pair).", seq_id))
            else:
                hepn1_center, hepn2_center = pair
                
                update_cursor.execute('''
                    UPDATE variants SET 
                        hepn1_start = ?, hepn1_end = ?, 
                        hepn2_start = ?, hepn2_end = ?, 
                        status = 'success', reason = 'HEPN anchored'
                    WHERE sequence_id = ?
                ''', (
                    max(0, hepn1_center - 30), hepn1_center + 80,
                    max(hepn1_center + 80, hepn2_center - 30), hepn2_center + 80,
                    seq_id
                ))
            processed += 1
            
        conn.commit()
        log.info(f"  Processed {processed} sequences.")

def generate_protenix_jsons(conn):
    """Builds Protenix JSONs from successful DB entries by chunking."""
    os.makedirs(JSON_OUT_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(METADATA_OUT_FILE), exist_ok=True)
    
    cursor = conn.cursor()
    cursor.execute('''
        SELECT sequence_id, sequence, crrna_repeat, hepn1_start, hepn1_end, hepn2_start, hepn2_end 
        FROM variants 
        WHERE status = 'success' AND crrna_repeat IS NOT NULL
    ''')
    
    domain_metadata_export = {}
    generated_count = 0
    
    print("Generating Protenix JSONs...")
    
    while True:
        batch = cursor.fetchmany(BATCH_SIZE)
        if not batch:
            break
            
        for row in batch:
            seq_id, protein_seq, dna_repeat, h1_start, h1_end, h2_start, h2_end = row
            
            # Convert DNA repeat to RNA
            dr_rna = dna_repeat.replace("T", "U").replace("t", "u")

            subtype = "unknown"
            if _HAS_CRRNA_HELPERS:
                spacer = _get_spacer(subtype)
                target_rna = _get_target(spacer)
                crrna_seq = _assemble_crrna(dr_rna, spacer, subtype)
            else:
                crrna_seq = dr_rna + DUMMY_SPACER_RNA
                target_rna = DUMMY_TARGET_RNA
            
            # Protenix expects proteinChain/rnaSequence (not protein/rna)
            protenix_payload = [
                {
                    "name": seq_id,
                    "sequences": [
                        {"proteinChain": {"sequence": protein_seq, "count": 1}},
                        {"rnaSequence": {"sequence": crrna_seq, "count": 1}},
                        {"rnaSequence": {"sequence": target_rna, "count": 1}}
                    ]
                }
            ]
            
            # Write JSON payload
            out_filename = os.path.join(JSON_OUT_DIR, f"{seq_id}.json")
            with open(out_filename, 'w') as out_f:
                json.dump(protenix_payload, out_f, indent=2)
                
            # Save exact HEPN windows to metadata
            domain_metadata_export[seq_id] = {
                "sequence_length": len(protein_seq),
                "domains": {
                    "HEPN1": {"start": h1_start, "end": h1_end},
                    "HEPN2": {"start": h2_start, "end": h2_end}
                },
                "crRNA_repeat_used": dr_rna,
                "subtype": "unknown",
            }
            generated_count += 1

    # Dump a flat JSON file of the HEPN metadata for downstream Phase 3
    with open(METADATA_OUT_FILE, 'w') as meta_f:
        json.dump(domain_metadata_export, meta_f, indent=2)

    return generated_count

_RUST_INGEST_BIN = shutil.which("cascade_ingest")
if not _RUST_INGEST_BIN:
    _candidate = os.path.join(os.path.dirname(__file__), "..", "rust", "target", "release", "cascade_ingest")
    if os.name == "nt":
        _candidate += ".exe"
    if os.path.isfile(_candidate):
        _RUST_INGEST_BIN = _candidate


def _try_rust_ingest():
    """Attempt to run the Rust-accelerated ingest pipeline. Returns True on success."""
    if not _RUST_INGEST_BIN:
        return False
    try:
        result = subprocess.run(
            [_RUST_INGEST_BIN,
             "--data-dir", DATA_DIR,
             "--db-file", DB_FILE,
             "--json-out-dir", JSON_OUT_DIR,
             "--metadata-out-file", METADATA_OUT_FILE],
            capture_output=False, text=True, timeout=600,
        )
        return result.returncode == 0
    except Exception as e:
        log.warning(f"Rust accelerator failed ({e}), falling back to Python.")
        return False


def main():
    if _try_rust_ingest():
        log.info("Rust-accelerated ingest completed successfully.")
        return

    log.info("Initializing SwitchBlade-Cas13 Local SQLite DB...")
    conn = init_db()
    
    log.info(f"Scanning {DATA_DIR} for fasta and csv files...")
    total_loaded = load_files_to_db(conn)
    log.info(f"Parsed {total_loaded} sequence boundaries into SQLite.")
    
    identify_hepn_domains(conn)
    
    generated = generate_protenix_jsons(conn)
    print(f"Successfully generated {generated} Protenix JSONs.")
    print(f"HEPN domain metadata saved to: {METADATA_OUT_FILE}")
    
    conn.close()

if __name__ == "__main__":
    main()