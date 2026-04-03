use clap::{Parser, Subcommand};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};
use std::process;
use walkdir::WalkDir;

#[derive(Parser)]
#[command(
    name = "cascade_structscore",
    about = "Structure file parsing and scoring for CASCADE pipeline"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Calculate Euclidean distance between two HEPN catalytic histidine CA atoms
    HepnDistance {
        #[arg(long)]
        structure: PathBuf,
        #[arg(long)]
        h1_idx: isize,
        #[arg(long)]
        h2_idx: isize,
        #[arg(long, default_value = "A")]
        chain: String,
    },
    /// Extract prediction scores from a Protenix summary JSON file
    ExtractScores {
        #[arg(long)]
        summary: PathBuf,
    },
    /// Recursively find structure files (.cif / .pdb) in a directory
    FindStructures {
        #[arg(long)]
        dir: PathBuf,
    },
}

fn load_structure(path: &Path) -> Result<pdbtbx::PDB, String> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_lowercase())
        .unwrap_or_default();

    let path_str = path.to_str().ok_or("Invalid path encoding")?;

    let result = match ext.as_str() {
        "cif" => pdbtbx::open_mmcif(path_str, pdbtbx::StrictnessLevel::Loose),
        "pdb" => pdbtbx::open_pdb(path_str, pdbtbx::StrictnessLevel::Loose),
        _ => return Err(format!("Unsupported file extension: .{ext}")),
    };

    result
        .map(|(pdb, _warnings)| pdb)
        .map_err(|errors| format!("{errors:?}"))
}

fn find_ca_coords(
    chain: &pdbtbx::Chain,
    residue_idx: isize,
    structure_path: &Path,
) -> Result<(f64, f64, f64), String> {
    let residue = chain
        .residues()
        .find(|r| r.serial_number() == residue_idx)
        .ok_or_else(|| {
            format!(
                "Could not find required residue or chain in structure {}. Error: \
                 residue {residue_idx} not found in chain {}",
                structure_path.display(),
                chain.id()
            )
        })?;

    let ca = residue
        .atoms()
        .find(|a| a.name().trim() == "CA")
        .ok_or_else(|| {
            format!(
                "Could not find required residue or chain in structure {}. Error: \
                 CA atom not found in residue {residue_idx}",
                structure_path.display()
            )
        })?;

    Ok(ca.pos())
}

fn cmd_hepn_distance(structure: &Path, h1_idx: isize, h2_idx: isize, chain_id: &str) {
    if !structure.exists() {
        eprintln!("Structure file not found: {}", structure.display());
        process::exit(1);
    }

    let pdb = match load_structure(structure) {
        Ok(pdb) => pdb,
        Err(e) => {
            eprintln!(
                "Could not find required residue or chain in structure {}. Error: {e}",
                structure.display()
            );
            process::exit(1);
        }
    };

    let model = match pdb.models().next() {
        Some(m) => m,
        None => {
            eprintln!(
                "Could not find required residue or chain in structure {}. Error: no models in file",
                structure.display()
            );
            process::exit(1);
        }
    };

    let chain = match model.chains().find(|c| c.id() == chain_id) {
        Some(c) => c,
        None => {
            eprintln!(
                "Could not find required residue or chain in structure {}. Error: chain '{chain_id}' not found",
                structure.display()
            );
            process::exit(1);
        }
    };

    let (x1, y1, z1) = match find_ca_coords(chain, h1_idx, structure) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("{e}");
            process::exit(1);
        }
    };

    let (x2, y2, z2) = match find_ca_coords(chain, h2_idx, structure) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("{e}");
            process::exit(1);
        }
    };

    let distance = ((x1 - x2).powi(2) + (y1 - y2).powi(2) + (z1 - z2).powi(2)).sqrt();
    println!("{:.3}", distance);
}

fn cmd_extract_scores(summary: &Path) {
    let contents = match std::fs::read_to_string(summary) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Failed to read summary file {}: {e}", summary.display());
            process::exit(1);
        }
    };

    let data: Value = match serde_json::from_str(&contents) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("Failed to parse JSON from {}: {e}", summary.display());
            process::exit(1);
        }
    };

    let iptm = data.get("iptm").and_then(Value::as_f64).unwrap_or(0.0);
    let ptm = data.get("ptm").and_then(Value::as_f64).unwrap_or(0.0);
    let ranking_score = data
        .get("ranking_score")
        .and_then(Value::as_f64)
        .unwrap_or(0.0);
    let af2_ig = data
        .get("af2_ig")
        .and_then(Value::as_f64)
        .or_else(|| data.get("af2_ig_score").and_then(Value::as_f64))
        .unwrap_or(0.0);

    let output = json!({
        "iptm": iptm,
        "ptm": ptm,
        "ranking_score": ranking_score,
        "af2_ig": af2_ig,
    });

    println!("{}", serde_json::to_string(&output).unwrap());
}

fn cmd_find_structures(dir: &Path) {
    if !dir.exists() || !dir.is_dir() {
        println!("[]");
        return;
    }

    let mut cif_files: Vec<String> = Vec::new();
    let mut pdb_files: Vec<String> = Vec::new();

    for entry in WalkDir::new(dir).into_iter().filter_map(Result::ok) {
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
            match ext.to_lowercase().as_str() {
                "cif" => cif_files.push(path.to_string_lossy().into_owned()),
                "pdb" => pdb_files.push(path.to_string_lossy().into_owned()),
                _ => {}
            }
        }
    }

    cif_files.sort();
    pdb_files.sort();

    let result = if !cif_files.is_empty() {
        &cif_files
    } else {
        &pdb_files
    };
    println!("{}", serde_json::to_string(result).unwrap());
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::HepnDistance {
            structure,
            h1_idx,
            h2_idx,
            chain,
        } => cmd_hepn_distance(structure, *h1_idx, *h2_idx, chain),
        Commands::ExtractScores { summary } => cmd_extract_scores(summary),
        Commands::FindStructures { dir } => cmd_find_structures(dir),
    }
}
