use clap::{Parser, Subcommand};
use csv::{ReaderBuilder, WriterBuilder};
use regex::Regex;
use serde_json::Value;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process;

#[derive(Parser)]
#[command(name = "cascade_sequtils", about = "Sequence analysis utilities for CASCADE")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Extract mutations between baseline and variant sequences
    ExtractMutations {
        #[arg(long)]
        baseline: PathBuf,
        #[arg(long)]
        variant: PathBuf,
    },
    /// Find catalytic histidine indices in a FASTA sequence
    FindHistidines {
        #[arg(long)]
        fasta: PathBuf,
    },
    /// Validate CRISPR repeat k-mers with optional structure filtering
    ValidateRepeats {
        #[arg(long)]
        data_dir: PathBuf,
        #[arg(long)]
        output_dir: PathBuf,
    },
}

fn main() {
    let cli = Cli::parse();
    let result = match cli.command {
        Commands::ExtractMutations { baseline, variant } => extract_mutations(&baseline, &variant),
        Commands::FindHistidines { fasta } => find_histidines(&fasta),
        Commands::ValidateRepeats { data_dir, output_dir } => {
            validate_repeats(&data_dir, &output_dir)
        }
    };
    if let Err(e) = result {
        eprintln!("Error: {e}");
        process::exit(1);
    }
}

// ---------------------------------------------------------------------------
// extract-mutations
// ---------------------------------------------------------------------------

fn read_baseline_sequence(path: &PathBuf) -> Result<String, String> {
    let content = fs::read_to_string(path)
        .map_err(|e| format!("Cannot read baseline file {}: {e}", path.display()))?;

    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    if ext.eq_ignore_ascii_case("json") {
        let parsed: Value =
            serde_json::from_str(&content).map_err(|e| format!("Invalid JSON: {e}"))?;
        // JSON format: [{name: ..., sequences: [{proteinChain: {sequence: ...}}, ...]}]
        let seq_obj = parsed
            .get(0)
            .and_then(|d| d.get("sequences"))
            .and_then(|s| s.get(0))
            .ok_or("JSON missing [0].sequences[0]")?;
        let seq = seq_obj
            .get("proteinChain")
            .or_else(|| seq_obj.get("protein"))
            .and_then(|p| p.get("sequence"))
            .and_then(|s| s.as_str())
            .ok_or("JSON missing proteinChain/protein.sequence")?;
        Ok(seq.to_string())
    } else {
        Ok(read_fasta_no_upper(&content))
    }
}

fn read_fasta_no_upper(content: &str) -> String {
    content
        .lines()
        .filter(|l| !l.starts_with('>'))
        .map(|l| l.trim())
        .collect::<Vec<_>>()
        .join("")
}

fn read_fasta_all(content: &str) -> String {
    content
        .lines()
        .filter(|l| !l.starts_with('>'))
        .map(|l| l.trim())
        .collect::<Vec<_>>()
        .join("")
        .to_uppercase()
}

fn read_fasta_first(content: &str) -> String {
    let mut seq_parts: Vec<&str> = Vec::new();
    let mut in_first = false;
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with('>') {
            if in_first {
                break;
            }
            in_first = true;
            continue;
        }
        if in_first || !trimmed.is_empty() {
            in_first = true;
            seq_parts.push(trimmed);
        }
    }
    seq_parts.join("").to_uppercase()
}

fn extract_mutations(baseline_path: &PathBuf, variant_path: &PathBuf) -> Result<(), String> {
    let baseline = read_baseline_sequence(baseline_path)?;
    let variant_content = fs::read_to_string(variant_path)
        .map_err(|e| format!("Cannot read variant file {}: {e}", variant_path.display()))?;
    let variant = read_fasta_no_upper(&variant_content);

    let max_len = baseline.len().max(variant.len());
    let baseline_bytes: Vec<u8> = {
        let mut v = baseline.into_bytes();
        v.resize(max_len, b'-');
        v
    };
    let variant_bytes: Vec<u8> = {
        let mut v = variant.into_bytes();
        v.resize(max_len, b'-');
        v
    };

    let mut mutations: Vec<String> = Vec::new();
    for i in 0..max_len {
        let b = baseline_bytes[i];
        let v = variant_bytes[i];
        if b != v {
            let pos = i + 1;
            if v == b'-' {
                mutations.push(format!("{pos}_del"));
            } else if b == b'-' {
                mutations.push(format!("{pos}_{}_ins", v as char));
            } else {
                mutations.push(format!("{pos}_{}", v as char));
            }
        }
    }

    let out = serde_json::to_string(&mutations).unwrap();
    println!("{out}");
    Ok(())
}

// ---------------------------------------------------------------------------
// find-histidines
// ---------------------------------------------------------------------------

fn find_histidines(fasta_path: &PathBuf) -> Result<(), String> {
    let content = fs::read_to_string(fasta_path)
        .map_err(|e| format!("Cannot read FASTA {}: {e}", fasta_path.display()))?;
    let seq = read_fasta_first(&content);

    let re = Regex::new(r"R.{3,6}H").unwrap();
    let matches: Vec<_> = re.find_iter(&seq).collect();

    if matches.len() < 2 {
        println!("{{\"h1\":null,\"h2\":null}}");
    } else {
        let h1 = matches[0].end() as i64;
        let h2 = matches[matches.len() - 1].end() as i64;
        println!("{{\"h1\":{h1},\"h2\":{h2}}}");
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// validate-repeats
// ---------------------------------------------------------------------------

const PREFERRED_MIN: usize = 28;
const PREFERRED_MAX: usize = 36;
const FALLBACK_MIN: usize = 23;
const FALLBACK_MAX: usize = 55;

fn is_simple_tandem_repeat(seq: &str, unit_max: usize, min_repeats: usize) -> bool {
    let upper = seq.to_uppercase();
    let bytes = upper.as_bytes();
    let slen = bytes.len();
    for unit_len in 2..=unit_max {
        let window = unit_len * min_repeats;
        if slen < window {
            continue;
        }
        for i in 0..=(slen - window) {
            let unit = &bytes[i..i + unit_len];
            let mut repeats: usize = 1;
            let mut pos = i + unit_len;
            while pos + unit_len <= slen && &bytes[pos..pos + unit_len] == unit {
                repeats += 1;
                pos += unit_len;
            }
            if repeats >= min_repeats
                && (repeats * unit_len) as f64 >= slen as f64 * 0.6
            {
                return true;
            }
        }
    }
    false
}

fn is_low_complexity(seq: &str, max_single_frac: f64) -> bool {
    if seq.is_empty() {
        return true;
    }
    let upper = seq.to_uppercase();
    let mut counts: HashMap<char, usize> = HashMap::new();
    for c in upper.chars() {
        *counts.entry(c).or_insert(0) += 1;
    }
    let max_count = counts.values().copied().max().unwrap_or(0);
    (max_count as f64 / seq.len() as f64) > max_single_frac
}

fn run_rnafold(sequence: &str) -> Option<(String, f64)> {
    let child = process::Command::new("RNAfold")
        .arg("--noPS")
        .stdin(process::Stdio::piped())
        .stdout(process::Stdio::piped())
        .stderr(process::Stdio::null())
        .spawn();

    let mut child = match child {
        Ok(c) => c,
        Err(_) => return None,
    };

    if let Some(ref mut stdin) = child.stdin {
        let _ = stdin.write_all(sequence.as_bytes());
        let _ = stdin.write_all(b"\n");
    }
    drop(child.stdin.take());

    let output = match child.wait_with_output() {
        Ok(o) => o,
        Err(_) => return None,
    };
    if !output.status.success() {
        return None;
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let last_line = stdout.lines().last()?;

    let mfe_re = Regex::new(r"\s*\([\s]*(-?\d+\.?\d*)\)\s*$").unwrap();
    let caps = mfe_re.captures(last_line)?;
    let mfe: f64 = caps.get(1)?.as_str().parse().ok()?;
    let mat = mfe_re.find(last_line)?;
    let structure = last_line[..mat.start()].trim().to_string();
    Some((structure, mfe))
}

fn has_plausible_stemloop(structure: &str, mfe: f64) -> bool {
    let paired: usize = structure.chars().filter(|&c| c == '(' || c == ')').count();
    if paired < 4 {
        return false;
    }
    let norm = mfe / structure.len() as f64;
    norm <= -0.05
}

fn validate_repeats(data_dir: &PathBuf, output_dir: &PathBuf) -> Result<(), String> {
    fs::create_dir_all(output_dir)
        .map_err(|e| format!("Cannot create output dir {}: {e}", output_dir.display()))?;

    let pattern = format!("{}/*_metadata.csv", data_dir.display());
    let mut paths: Vec<PathBuf> = glob::glob(&pattern)
        .map_err(|e| format!("Glob error: {e}"))?
        .filter_map(|r| r.ok())
        .collect();
    paths.sort();

    #[derive(Default)]
    struct Row {
        sequence_id: String,
        original_first_kmer: String,
        chosen_repeat: String,
        chosen_length: usize,
        selection_reason: String,
        structure: String,
        mfe_kcal_mol: String,
        structure_ok: bool,
    }

    let mut rows: Vec<Row> = Vec::new();

    for csv_path in &paths {
        let mut rdr = ReaderBuilder::new()
            .has_headers(true)
            .from_path(csv_path)
            .map_err(|e| format!("Cannot read CSV {}: {e}", csv_path.display()))?;

        let headers = rdr
            .headers()
            .map_err(|e| format!("Bad CSV headers: {e}"))?
            .clone();

        let seq_id_idx = headers.iter().position(|h| h == "sequence_id");
        let repeat_idx = headers.iter().position(|h| h == "repeat_domains");

        let seq_id_idx =
            seq_id_idx.ok_or_else(|| format!("Missing sequence_id column in {}", csv_path.display()))?;
        let repeat_idx =
            repeat_idx.ok_or_else(|| format!("Missing repeat_domains column in {}", csv_path.display()))?;

        for result in rdr.records() {
            let record = result.map_err(|e| format!("CSV read error: {e}"))?;
            let sequence_id = record.get(seq_id_idx).unwrap_or("").to_string();
            let repeat_domains_raw = record.get(repeat_idx).unwrap_or("").to_string();

            let raw_kmers: Vec<String> = repeat_domains_raw
                .split('|')
                .map(|s| s.trim().to_string())
                .collect();

            let original_first_kmer = raw_kmers
                .first()
                .map(|s| {
                    let chars: Vec<char> = s.chars().collect();
                    chars.iter().take(50).collect::<String>()
                })
                .unwrap_or_default();

            let candidates: Vec<String> = raw_kmers
                .iter()
                .map(|k| k.replace('T', "U").replace('t', "u"))
                .filter(|k| {
                    !k.is_empty()
                        && k.len() >= FALLBACK_MIN
                        && k.len() <= FALLBACK_MAX
                        && !is_simple_tandem_repeat(k, 6, 3)
                        && !is_low_complexity(k, 0.5)
                })
                .collect();

            let preferred: Vec<&String> = candidates
                .iter()
                .filter(|k| k.len() >= PREFERRED_MIN && k.len() <= PREFERRED_MAX)
                .collect();

            let (chosen, reason) = if !preferred.is_empty() {
                let best = preferred.iter().max_by_key(|k| k.len()).unwrap();
                (Some((*best).clone()), format!("preferred_len_{}nt", best.len()))
            } else if !candidates.is_empty() {
                let best = candidates.iter().max_by_key(|k| k.len()).unwrap();
                (Some(best.clone()), format!("fallback_longest_{}nt", best.len()))
            } else {
                (None, "no_valid_candidates".to_string())
            };

            let mut row = Row {
                sequence_id,
                original_first_kmer,
                chosen_repeat: String::new(),
                chosen_length: 0,
                selection_reason: reason,
                structure: String::new(),
                mfe_kcal_mol: String::new(),
                structure_ok: false,
            };

            if let Some(ref seq) = chosen {
                row.chosen_repeat = seq.clone();
                row.chosen_length = seq.len();

                if let Some((structure, mfe)) = run_rnafold(seq) {
                    row.structure_ok = has_plausible_stemloop(&structure, mfe);
                    row.structure = structure;
                    row.mfe_kcal_mol = format!("{mfe:.2}");
                }
            }

            rows.push(row);
        }
    }

    let report_path = output_dir.join("repeat_validation_report.csv");
    let ids_path = output_dir.join("validated_baseline_ids.txt");

    {
        let mut wtr = WriterBuilder::new()
            .from_path(&report_path)
            .map_err(|e| format!("Cannot write report CSV: {e}"))?;

        wtr.write_record([
            "sequence_id",
            "original_first_kmer",
            "chosen_repeat",
            "chosen_length",
            "selection_reason",
            "structure",
            "mfe_kcal_mol",
            "structure_ok",
        ])
        .map_err(|e| format!("CSV write error: {e}"))?;

        for r in &rows {
            wtr.write_record([
                &r.sequence_id,
                &r.original_first_kmer,
                &r.chosen_repeat,
                &r.chosen_length.to_string(),
                &r.selection_reason,
                &r.structure,
                &r.mfe_kcal_mol,
                &if r.structure_ok {
                    "true".to_string()
                } else {
                    "false".to_string()
                },
            ])
            .map_err(|e| format!("CSV write error: {e}"))?;
        }
        wtr.flush().map_err(|e| format!("CSV flush error: {e}"))?;
    }

    let validated_ids: Vec<&str> = rows
        .iter()
        .filter(|r| r.structure_ok)
        .map(|r| r.sequence_id.as_str())
        .collect();

    {
        let mut f = fs::File::create(&ids_path)
            .map_err(|e| format!("Cannot write IDs file: {e}"))?;
        for id in &validated_ids {
            writeln!(f, "{id}").map_err(|e| format!("Write error: {e}"))?;
        }
    }

    let total = rows.len();
    let passed = validated_ids.len();
    let has_chosen = rows.iter().any(|r| !r.chosen_repeat.is_empty());

    println!("Wrote {total} rows to {}", report_path.display());
    println!("Passed structure filter: {passed}/{total}");
    if passed == 0 && has_chosen {
        println!("  NOTE: Install ViennaRNA for structure filtering: pip install ViennaRNA");
    }
    println!("Validated IDs written to {}", ids_path.display());
    if !validated_ids.is_empty() {
        let show: Vec<&str> = validated_ids.iter().take(5).copied().collect();
        println!("  First 5: {}", show.join(", "));
    }

    Ok(())
}
