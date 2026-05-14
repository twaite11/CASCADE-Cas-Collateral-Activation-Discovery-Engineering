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
    /// Calculate Euclidean distance between two HEPN catalytic histidine side-chain atoms.
    ///
    /// By default the distance is measured between the imidazole NE2 nitrogens
    /// (the catalytic atoms that coordinate the scissile-phosphate water in the
    /// composite HEPN active site).  The CASCADE OFF >= 18 A and ON <= 12 A gates
    /// are defined on this NE2-NE2 distance (see scripts/utils/pdb_kinematics.py).
    /// If NE2 is missing (mutated to a non-His residue, or model lacks side-chain
    /// atoms), the tool falls back through the comma-separated --atom-fallback
    /// list (default: ND1, CG, CA) so a number is still returned but with a
    /// well-defined provenance.
    HepnDistance {
        #[arg(long)]
        structure: PathBuf,
        #[arg(long)]
        h1_idx: isize,
        #[arg(long)]
        h2_idx: isize,
        #[arg(long, default_value = "A")]
        chain: String,
        #[arg(long, default_value = "NE2")]
        atom: String,
        #[arg(long, default_value = "ND1,CG,CA")]
        atom_fallback: String,
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

fn find_atom_coords(
    chain: &pdbtbx::Chain,
    residue_idx: isize,
    atom_chain: &[&str],
    structure_path: &Path,
) -> Result<((f64, f64, f64), String), String> {
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

    for name in atom_chain {
        if let Some(atom) = residue.atoms().find(|a| a.name().trim() == *name) {
            return Ok((atom.pos(), (*name).to_string()));
        }
    }

    Err(format!(
        "Could not find required residue or chain in structure {}. Error: \
         none of atoms [{}] found in residue {residue_idx}",
        structure_path.display(),
        atom_chain.join(",")
    ))
}

fn cmd_hepn_distance(
    structure: &Path,
    h1_idx: isize,
    h2_idx: isize,
    chain_id: &str,
    atom_primary: &str,
    atom_fallback_csv: &str,
) {
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

    // Build the atom-name preference list: primary first, then comma-separated
    // fallbacks (defaults NE2 -> ND1 -> CG -> CA).  NE2 is the catalytic atom
    // in the Cas13 composite HEPN active site.
    let mut atom_chain: Vec<&str> = Vec::with_capacity(8);
    atom_chain.push(atom_primary);
    for name in atom_fallback_csv.split(',') {
        let trimmed = name.trim();
        if !trimmed.is_empty() && trimmed != atom_primary {
            atom_chain.push(trimmed);
        }
    }

    let (p1, a1) = match find_atom_coords(chain, h1_idx, &atom_chain, structure) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("{e}");
            process::exit(1);
        }
    };
    let (p2, a2) = match find_atom_coords(chain, h2_idx, &atom_chain, structure) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("{e}");
            process::exit(1);
        }
    };

    let (x1, y1, z1) = p1;
    let (x2, y2, z2) = p2;
    let distance = ((x1 - x2).powi(2) + (y1 - y2).powi(2) + (z1 - z2).powi(2)).sqrt();

    // Print distance only on stdout (preserves the existing CLI contract: a
    // single floating-point number).  Provenance of which atoms were measured
    // goes to stderr so consumers parsing stdout aren't affected.
    if a1 != atom_primary || a2 != atom_primary {
        eprintln!(
            "[cascade_structscore] hepn-distance fell back: h1={a1} h2={a2} (primary={atom_primary})"
        );
    }
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
    let mut af2_ig = data
        .get("af2_ig")
        .and_then(Value::as_f64)
        .or_else(|| data.get("af2_ig_score").and_then(Value::as_f64))
        .unwrap_or(0.0);

    // Derive AF2-IG from chain_pair_iptm off-diagonal mean when explicit key is absent
    if af2_ig == 0.0 {
        if let Some(cp_matrix) = data.get("chain_pair_iptm").and_then(Value::as_array) {
            let n = cp_matrix.len();
            let mut sum = 0.0_f64;
            let mut count = 0_u64;
            for (i, row) in cp_matrix.iter().enumerate() {
                if let Some(row_arr) = row.as_array() {
                    for (j, val) in row_arr.iter().enumerate() {
                        if i != j {
                            if let Some(v) = val.as_f64() {
                                sum += v;
                                count += 1;
                            }
                        }
                    }
                }
            }
            if count > 0 {
                af2_ig = sum / count as f64;
            }
        }
        // Fallback: weighted proxy from global scores
        if af2_ig == 0.0 && (iptm > 0.0 || ptm > 0.0) {
            af2_ig = 0.8 * iptm + 0.2 * ptm;
        }
    }

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
            atom,
            atom_fallback,
        } => cmd_hepn_distance(structure, *h1_idx, *h2_idx, chain, atom, atom_fallback),
        Commands::ExtractScores { summary } => cmd_extract_scores(summary),
        Commands::FindStructures { dir } => cmd_find_structures(dir),
    }
}
