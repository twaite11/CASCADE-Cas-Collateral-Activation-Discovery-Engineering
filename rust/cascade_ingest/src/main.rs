use std::collections::BTreeMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process;
use std::time::SystemTime;

use clap::Parser;
use regex::Regex;
use rusqlite::{params, Connection};
use serde_json::{json, Map, Value};

#[derive(Parser)]
#[command(name = "cascade_ingest", about = "Fast FASTA/CSV ingest and HEPN scanning for CASCADE")]
struct Cli {
    #[arg(long)]
    data_dir: String,
    #[arg(long)]
    db_file: String,
    #[arg(long)]
    json_out_dir: String,
    #[arg(long)]
    metadata_out_file: String,
}

fn timestamp() -> String {
    let now = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default();
    let total_secs = now.as_secs();
    let secs_of_day = total_secs % 86400;
    let h = secs_of_day / 3600;
    let m = (secs_of_day % 3600) / 60;
    let s = secs_of_day % 60;
    format!("{:02}:{:02}:{:02}", h, m, s)
}

fn log_info(msg: &str) {
    eprintln!("[{}] INFO {}", timestamp(), msg);
}

fn ensure_parent_dir(path: &str) -> Result<(), String> {
    if let Some(parent) = Path::new(path).parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .map_err(|e| format!("Failed to create directory {}: {}", parent.display(), e))?;
        }
    }
    Ok(())
}

fn sorted_glob(pattern: &str) -> Vec<std::path::PathBuf> {
    let mut paths: Vec<_> = glob::glob(pattern)
        .unwrap_or_else(|e| {
            eprintln!("Bad glob pattern '{}': {}", pattern, e);
            process::exit(1);
        })
        .filter_map(Result::ok)
        .collect();
    paths.sort();
    paths
}

fn run() -> Result<(), String> {
    let cli = Cli::parse();

    // Step 1: Init DB
    log_info("Initializing SwitchBlade-Cas13 Local SQLite DB...");
    ensure_parent_dir(&cli.db_file)?;
    let mut conn = Connection::open(&cli.db_file)
        .map_err(|e| format!("Failed to open database: {}", e))?;
    conn.execute_batch("PRAGMA journal_mode=WAL;")
        .map_err(|e| format!("Failed to set journal mode: {}", e))?;
    conn.execute(
        "CREATE TABLE IF NOT EXISTS variants (
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
        )",
        [],
    )
    .map_err(|e| format!("Failed to create table: {}", e))?;

    // Step 2: Parse FASTA files
    log_info("Streaming FASTA sequences to SQLite...");
    let fasta_pattern = format!("{}/*.fasta", cli.data_dir);
    let fasta_files = sorted_glob(&fasta_pattern);

    for fasta_path in &fasta_files {
        let file = fs::File::open(fasta_path)
            .map_err(|e| format!("Failed to open {}: {}", fasta_path.display(), e))?;
        let reader = BufReader::new(file);

        let tx = conn.transaction().map_err(|e| format!("Transaction error: {}", e))?;
        {
            let mut current_id: Option<String> = None;
            let mut current_seq = String::new();
            let mut stmt = tx
                .prepare(
                    "INSERT INTO variants (sequence_id, sequence) VALUES (?1, ?2)
                     ON CONFLICT(sequence_id) DO UPDATE SET sequence=excluded.sequence",
                )
                .map_err(|e| format!("Failed to prepare FASTA upsert: {}", e))?;

            for line in reader.lines() {
                let line = line.map_err(|e| format!("Read error: {}", e))?;
                if let Some(header) = line.strip_prefix('>') {
                    if let Some(ref id) = current_id {
                        stmt.execute(params![id, current_seq])
                            .map_err(|e| format!("Failed to upsert sequence: {}", e))?;
                    }
                    current_id = Some(header.to_string());
                    current_seq.clear();
                } else {
                    current_seq.push_str(line.trim());
                }
            }
            if let Some(ref id) = current_id {
                stmt.execute(params![id, current_seq])
                    .map_err(|e| format!("Failed to upsert sequence: {}", e))?;
            }
        }
        tx.commit().map_err(|e| format!("Commit error: {}", e))?;
    }

    // Step 3: Parse CSV files
    log_info("Streaming CSV metadata to SQLite...");
    let csv_pattern = format!("{}/*.csv", cli.data_dir);
    let csv_files = sorted_glob(&csv_pattern);

    for csv_path in &csv_files {
        let mut rdr = csv::Reader::from_path(csv_path)
            .map_err(|e| format!("Failed to open CSV {}: {}", csv_path.display(), e))?;

        let headers = rdr
            .headers()
            .map_err(|e| format!("CSV header error: {}", e))?
            .clone();

        let tx = conn.transaction().map_err(|e| format!("Transaction error: {}", e))?;
        {
            let mut stmt = tx
                .prepare(
                    "UPDATE variants SET crrna_repeat = ?1, sra_accession = ?2, score = ?3
                     WHERE sequence_id = ?4",
                )
                .map_err(|e| format!("Failed to prepare CSV update: {}", e))?;

            for result in rdr.records() {
                let record = result.map_err(|e| format!("CSV parse error: {}", e))?;

                let idx = |name: &str| -> Option<usize> {
                    headers.iter().position(|h| h == name)
                };

                let seq_id = idx("sequence_id")
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .to_string();
                let repeat_raw = idx("repeat_domains")
                    .and_then(|i| record.get(i))
                    .unwrap_or("");
                let repeat = repeat_raw.split('|').next().unwrap_or("").trim().to_string();
                let sra = idx("sra_accession")
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .to_string();
                let score: f64 = idx("score")
                    .and_then(|i| record.get(i))
                    .unwrap_or("0.0")
                    .parse()
                    .unwrap_or(0.0);

                stmt.execute(params![repeat, sra, score, seq_id])
                    .map_err(|e| format!("Failed to update CSV metadata: {}", e))?;
            }
        }
        tx.commit().map_err(|e| format!("Commit error: {}", e))?;
    }

    // Step 4: HEPN domain scanning
    log_info("Scanning sequences for HEPN domains (may take a minute for large DBs)...");
    let hepn_re = Regex::new(r"R.{3,6}H").unwrap();

    let rows: Vec<(String, String)> = {
        let mut stmt = conn
            .prepare("SELECT sequence_id, sequence FROM variants WHERE sequence IS NOT NULL")
            .map_err(|e| format!("Failed to prepare HEPN select: {}", e))?;
        stmt.query_map([], |row| {
            Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
        })
        .map_err(|e| format!("Failed to query sequences: {}", e))?
        .filter_map(Result::ok)
        .collect()
    };

    let mut processed = 0usize;
    for chunk in rows.chunks(1000) {
        let tx = conn.transaction().map_err(|e| format!("Transaction error: {}", e))?;
        for (seq_id, sequence) in chunk {
            let matches: Vec<_> = hepn_re.find_iter(sequence).collect();
            if matches.len() < 2 {
                let reason = format!("Only {} HEPN motifs found.", matches.len());
                tx.execute(
                    "UPDATE variants SET status = 'failed', reason = ?1 WHERE sequence_id = ?2",
                    params![reason, seq_id],
                )
                .map_err(|e| format!("Failed to update failed HEPN: {}", e))?;
            } else {
                let hepn1_center = matches[0].start() as i64;
                let hepn2_center = matches[matches.len() - 1].start() as i64;
                let h1_start = 0i64.max(hepn1_center - 30);
                let h1_end = hepn1_center + 80;
                let h2_start = (hepn1_center + 80).max(hepn2_center - 30);
                let h2_end = hepn2_center + 80;
                tx.execute(
                    "UPDATE variants SET hepn1_start=?1, hepn1_end=?2, hepn2_start=?3, hepn2_end=?4, \
                     status='success', reason='HEPN anchored' WHERE sequence_id=?5",
                    params![h1_start, h1_end, h2_start, h2_end, seq_id],
                )
                .map_err(|e| format!("Failed to update HEPN domains: {}", e))?;
            }
        }
        tx.commit().map_err(|e| format!("Commit error: {}", e))?;
        processed += chunk.len();
        eprintln!("  Processed {} sequences.", processed);
    }

    // Step 5: Generate Protenix JSONs
    println!("Generating Protenix JSONs...");
    fs::create_dir_all(&cli.json_out_dir)
        .map_err(|e| format!("Failed to create json output dir: {}", e))?;
    ensure_parent_dir(&cli.metadata_out_file)?;

    const DUMMY_SPACER_RNA: &str = "GUCGACUGACGUACGUACGUACGU";
    const DUMMY_TARGET_RNA: &str = "AAAAAACGUACGUACGUACGUCAGUCGACAAAAAA";

    let success_rows: Vec<(String, String, String, i64, i64, i64, i64)> = {
        let mut stmt = conn
            .prepare(
                "SELECT sequence_id, sequence, crrna_repeat, hepn1_start, hepn1_end, hepn2_start, hepn2_end \
                 FROM variants WHERE status = 'success' AND crrna_repeat IS NOT NULL",
            )
            .map_err(|e| format!("Failed to prepare Protenix select: {}", e))?;
        stmt.query_map([], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, String>(1)?,
                row.get::<_, String>(2)?,
                row.get::<_, i64>(3)?,
                row.get::<_, i64>(4)?,
                row.get::<_, i64>(5)?,
                row.get::<_, i64>(6)?,
            ))
        })
        .map_err(|e| format!("Failed to query success rows: {}", e))?
        .filter_map(Result::ok)
        .collect()
    };

    let mut metadata: BTreeMap<String, Value> = BTreeMap::new();
    let mut count = 0usize;

    for (seq_id, protein_seq, crrna_repeat, h1_start, h1_end, h2_start, h2_end) in &success_rows {
        let dr_rna = crrna_repeat.replace('T', "U").replace('t', "u");
        let crrna_seq = format!("{}{}", dr_rna, DUMMY_SPACER_RNA);

        let payload = json!([{
            "name": seq_id,
            "sequences": [
                {"proteinChain": {"sequence": protein_seq, "count": 1}},
                {"rnaSequence": {"sequence": crrna_seq, "count": 1}},
                {"rnaSequence": {"sequence": DUMMY_TARGET_RNA, "count": 1}}
            ]
        }]);

        let json_path = Path::new(&cli.json_out_dir).join(format!("{}.json", seq_id));
        let formatted = serde_json::to_string_pretty(&payload)
            .map_err(|e| format!("JSON serialize error: {}", e))?;
        fs::write(&json_path, format!("{}\n", formatted))
            .map_err(|e| format!("Failed to write {}: {}", json_path.display(), e))?;

        let mut domains = Map::new();
        domains.insert("HEPN1".to_string(), json!({"start": h1_start, "end": h1_end}));
        domains.insert("HEPN2".to_string(), json!({"start": h2_start, "end": h2_end}));

        let entry = json!({
            "sequence_length": protein_seq.len(),
            "domains": domains,
            "crRNA_repeat_used": dr_rna
        });
        metadata.insert(seq_id.clone(), entry);
        count += 1;
    }

    println!("Successfully generated {} Protenix JSONs.", count);

    let meta_json = serde_json::to_string_pretty(&metadata)
        .map_err(|e| format!("Metadata JSON error: {}", e))?;
    fs::write(&cli.metadata_out_file, format!("{}\n", meta_json))
        .map_err(|e| format!("Failed to write metadata: {}", e))?;

    println!("HEPN domain metadata saved to: {}", cli.metadata_out_file);

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}
