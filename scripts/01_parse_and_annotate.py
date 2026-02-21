import csv
import json
import os
import re
import glob
import sqlite3
import logging

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

# RNA Constants
DUMMY_SPACER_RNA = "GUCGACUGACGUACGUACGUACGU" # 24nt
DUMMY_TARGET_RNA = "AAAAAA" + "ACGUACGUACGUACGUCAGUCGAC" + "AAAAAA" # Simulates tumor fusion RNA
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
    for csv_file in csv_files:
        with open(csv_file, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                seq_id = row.get("sequence_id")
                repeat_str = row.get("repeat_domains", "").split("|")[0]
                sra_acc = row.get("sra_accession", "")
                score = float(row.get("score", 0.0))
                
                cursor.execute('''
                    UPDATE variants 
                    SET crrna_repeat = ?, sra_accession = ?, score = ?
                    WHERE sequence_id = ?
                ''', (repeat_str, sra_acc, score, seq_id))
    conn.commit()
                
    return items_processed

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
    
    motif = re.compile(r'R[AILMFVWY][A-Z]{3}H')
    processed = 0
    
    while True:
        # Fetch a chunk of records to keep memory usage low
        batch = cursor.fetchmany(BATCH_SIZE)
        if not batch:
            break
            
        for seq_id, sequence in batch:
            matches = list(motif.finditer(sequence))
            
            if len(matches) < 2:
                update_cursor.execute('''
                    UPDATE variants SET status = 'failed', reason = ? WHERE sequence_id = ?
                ''', (f"Only {len(matches)} HEPN motifs found.", seq_id))
            else:
                hepn1_center = matches[0].start()
                hepn2_center = matches[-1].start()
                
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
            crrna_seq = dr_rna + DUMMY_SPACER_RNA
            
            protenix_payload = [
                {
                    "name": seq_id,
                    "sequences": [
                        {"protein": {"id": "A", "sequence": protein_seq}},
                        {"rna": {"id": "B", "sequence": crrna_seq}},
                        {"rna": {"id": "C", "sequence": DUMMY_TARGET_RNA}}
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
                "crRNA_repeat_used": dr_rna
            }
            generated_count += 1

    # Dump a flat JSON file of the HEPN metadata for downstream Phase 3
    with open(METADATA_OUT_FILE, 'w') as meta_f:
        json.dump(domain_metadata_export, meta_f, indent=2)

    return generated_count

def main():
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