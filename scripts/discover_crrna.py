#!/usr/bin/env python3
"""
Iterative crRNA Direct Repeat (DR) Discovery Pipeline for novel Cas13 candidates.

Exhausts 5 tiers of increasingly sophisticated approaches:
  Tier 1: Proximity CRISPR detection near the Cas13 gene on source contigs
  Tier 2: Salvage/validate DRs from original mining repeat_domains
  Tier 3: Known Cas13 DR library screen via Protenix ipTM
  Tier 4: Broader metagenomic CRISPR array search
  Tier 5: Computational DR design (stem-loop generation + mutagenesis)

Each tier produces candidate DRs scored by Protenix mini ipTM.
Short-circuits for any enzyme achieving ipTM > threshold (default 0.5).

Usage:
  python discover_crrna.py --contigs /tmp/all_contigs.fasta
  python discover_crrna.py --contigs /tmp/all_contigs.fasta --skip-fold  # DR discovery only
  python discover_crrna.py --tier 3  # start from tier 3
"""
import argparse
import csv
import json
import logging
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
METADATA_JSON = PROJECT_ROOT / "metadata" / "variant_domain_metadata.json"
DB_FILE = PROJECT_ROOT / "metadata" / "cas13_variants.db"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
REPORT_CSV = OUTPUT_DIR / "crrna_discovery_report.csv"
JSONS_DIR = PROJECT_ROOT / "jsons"
STRUCTURE_CHECK_DIR = OUTPUT_DIR / "structure_check"

sys.path.insert(0, str(SCRIPT_DIR))
from fix_crrna_assignments import (
    find_crispr_arrays,
    is_trna_like,
    has_crispr_structure,
    _pick_best_dr_orientation,
    _salvage_from_old_kmers,
)
from utils.protenix_eval import (
    run_protenix_inference,
    assemble_crrna,
    get_spacer_for_subtype,
    get_target_for_spacer,
)
from utils.pdb_kinematics import extract_protenix_scores

IPTM_THRESHOLD = 0.5

CANDIDATES = {
    "3174394687_r3_1090aa_rev": {
        "contig": "3174394687",
        "subtype": "cas13a",
        "blastx_qstart": 291032,
        "blastx_qend": 290433,
    },
    "3174386807_r2_736aa_rev": {
        "contig": "3174386807",
        "subtype": "cas13a",
        "blastx_qstart": 64091,
        "blastx_qend": 63492,
    },
    "3174387268_f1_1023aa_fwd": {
        "contig": "3174387268",
        "subtype": "cas13b",
        "blastx_qstart": 155125,
        "blastx_qend": 154673,
    },
    "3174383916_f3_920aa_fwd": {
        "contig": "3174383916",
        "subtype": "cas13a",
        "blastx_qstart": 144673,
        "blastx_qend": 144416,
    },
}

# ---------------------------------------------------------------------------
# Known Cas13 DR library (Tier 3) — DNA sequences, converted to RNA at runtime
# ---------------------------------------------------------------------------
KNOWN_DR_LIBRARY = {
    "LshCas13a":  "GATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC",
    "LbuCas13a":  "GACCACCCCAAAAATGAAGGGGACTAAAAC",
    "LwaCas13a":  "GATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC",
    "LseCas13a":  "GATTTAGATTAACCCCTCAAAAGGGACTAAAAT",
    "HheCas13a":  "GATTTAGACCACCCCAAAAAATGAAGGGGACTAAAAC",
    "PbuCas13b":  "GTTTTGATAAACCATTAATAAGATTGATTTTAAACC",
    "BzCas13b":   "GTTTTGATAAACTTAATCAAGTCTATTTGAAACC",
    "PsmCas13b":  "GTTTTGTTATAAAACGTTTTGAAATTTCC",
    "Pin2Cas13b": "GTTGTTGAAATCCTCCCTTAGAGGGATTTAAC",
    "RfxCas13d":  "AACCCCCACCCCGCGGGGGGATTTTTTTAT",
    "EsCas13d":   "CAACCATATCCCCCTATCGCAGGGATTTTTAT",
    "AdmCas13d":  "CAACCATACCATTGTATGCAGGGATTTTTTTAT",
    "UrCas13d":   "AACCCCTACCCCGCGAGGGGGATTTTTTAT",
    "Cas13bt1":   "GTTTCGAAATTTCATTAAACTTGACC",
    "Cas13bt3":   "GTTCATTGATATTCGTTACGCTGATTTAAAC",
}


# =====================================================================
# Shared: Protenix ipTM scoring
# =====================================================================

def _load_protein_seq(candidate_id: str) -> str:
    """Load protein sequence from SQLite DB."""
    conn = sqlite3.connect(str(DB_FILE))
    cur = conn.cursor()
    cur.execute("SELECT sequence FROM variants WHERE sequence_id=?", (candidate_id,))
    row = cur.fetchone()
    conn.close()
    if not row:
        raise ValueError(f"Candidate {candidate_id} not found in {DB_FILE}")
    return row[0]


def score_dr_with_protenix(
    candidate_id: str,
    protein_seq: str,
    dr_rna: str,
    subtype: str,
    tag: str = "",
    skip_fold: bool = False,
) -> dict:
    """Score a DR candidate by running Protenix mini fold and extracting ipTM.

    Returns dict with keys: dr_rna, iptm, ptm, af2_ig, tag, json_path, summary_path
    """
    spacer = get_spacer_for_subtype(subtype)
    target_rna = get_target_for_spacer(spacer)
    crrna = assemble_crrna(dr_rna, spacer, subtype)

    run_tag = f"{candidate_id}_{tag}" if tag else candidate_id
    safe_tag = re.sub(r"[^\w\-.]", "_", run_tag)

    out_base = STRUCTURE_CHECK_DIR / "dr_screen"
    out_base.mkdir(parents=True, exist_ok=True)
    json_path = out_base / f"{safe_tag}.json"

    payload = [
        {
            "name": safe_tag,
            "sequences": [
                {"proteinChain": {"sequence": protein_seq, "count": 1}},
                {"rnaSequence": {"sequence": crrna, "count": 1}},
            ],
        }
    ]
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    result = {
        "candidate_id": candidate_id,
        "dr_rna": dr_rna,
        "dr_len": len(dr_rna),
        "tag": tag,
        "iptm": 0.0,
        "ptm": 0.0,
        "af2_ig": 0.0,
        "json_path": str(json_path),
        "summary_path": "",
    }

    if skip_fold:
        return result

    try:
        out_dir = str(out_base / f"{safe_tag}_results")
        struct_path, summary_path = run_protenix_inference(
            json_path=str(json_path),
            out_dir=out_dir,
            model_tier="mini",
        )
        result["summary_path"] = str(summary_path)
        scores = extract_protenix_scores(summary_path)
        result["iptm"] = scores.get("iptm", 0.0)
        result["ptm"] = scores.get("ptm", 0.0)
        result["af2_ig"] = scores.get("af2_ig", 0.0)
    except Exception as e:
        log.warning(f"Protenix scoring failed for {safe_tag}: {e}")

    return result


# =====================================================================
# Tier 1: Proximity-based CRISPR detection
# =====================================================================

def _load_contigs(contigs_path: str) -> dict:
    """Load contigs from FASTA into {id: seq} dict."""
    contigs = {}
    current_id = ""
    current_seq = []
    with open(contigs_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    contigs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        contigs[current_id] = "".join(current_seq)
    return contigs


def _reverse_complement_dna(seq: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return seq.upper().translate(comp)[::-1]


def _run_minced_on_seq(dna_seq: str, contig_id: str = "query") -> list:
    """Run minced on a single sequence, return list of DR candidates (RNA)."""
    minced_bin = shutil.which("minced")
    if not minced_bin:
        return []

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        tmp.write(f">{contig_id}\n{dna_seq}\n")
        tmp_path = tmp.name
    out_path = tmp_path + ".out"

    try:
        subprocess.run(
            [minced_bin, "-minNR", "2", "-minRL", "20", "-maxRL", "50", tmp_path, out_path],
            capture_output=True, text=True, timeout=120,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return []

    drs = []
    if os.path.exists(out_path):
        with open(out_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("POSITION") or line.startswith("-") or line.startswith("Repeats") or line.startswith("CRISPR") or line.startswith("Sequence") or line.startswith("Time"):
                    continue
                cols = line.split()
                if len(cols) < 2:
                    continue
                try:
                    int(cols[0])
                    repeat_candidate = cols[1]
                except (ValueError, IndexError):
                    continue
                if len(repeat_candidate) >= 20 and set(repeat_candidate.upper()) <= {"A", "C", "G", "T"}:
                    rna = repeat_candidate.upper().replace("T", "U")
                    drs.append(rna)

    for tmp_f in [tmp_path, out_path]:
        try:
            os.unlink(tmp_f)
        except OSError:
            pass

    return list(set(drs))


def tier1_proximity_detection(contigs: dict, skip_fold: bool = False) -> dict:
    """Tier 1: Search for CRISPR arrays in +-20kb window around each gene."""
    log.info("=" * 60)
    log.info("TIER 1: Proximity-based CRISPR detection")
    log.info("=" * 60)

    results = {}
    window = 20000

    for cid, info in CANDIDATES.items():
        contig_id = info["contig"]
        if contig_id not in contigs:
            log.warning(f"Contig {contig_id} not found for {cid}")
            results[cid] = []
            continue

        contig_seq = contigs[contig_id].upper()
        gene_center = (info["blastx_qstart"] + info["blastx_qend"]) // 2
        start = max(0, gene_center - window)
        end = min(len(contig_seq), gene_center + window)
        neighborhood = contig_seq[start:end]

        log.info(f"{cid}: searching {start}-{end} ({end - start}bp) around gene center {gene_center}")

        candidate_drs = []

        # Python detector on forward strand, relaxed to min_units=2
        arrays_fwd = find_crispr_arrays(neighborhood, min_units=2)
        for arr in arrays_fwd:
            dr = arr["consensus_repeat"]
            if not is_trna_like(dr):
                candidate_drs.append(("python_fwd", dr, arr.get("n_repeats", 0)))

        # Python detector on reverse complement
        rc_neighborhood = _reverse_complement_dna(neighborhood)
        arrays_rc = find_crispr_arrays(rc_neighborhood, min_units=2)
        for arr in arrays_rc:
            dr = arr["consensus_repeat"]
            if not is_trna_like(dr):
                candidate_drs.append(("python_rc", dr, arr.get("n_repeats", 0)))

        # MinCED on forward
        minced_drs = _run_minced_on_seq(neighborhood, f"{contig_id}_neighborhood")
        for dr in minced_drs:
            if not is_trna_like(dr):
                candidate_drs.append(("minced_fwd", dr, 0))

        # MinCED on reverse complement
        minced_rc_drs = _run_minced_on_seq(rc_neighborhood, f"{contig_id}_neighborhood_rc")
        for dr in minced_rc_drs:
            if not is_trna_like(dr):
                candidate_drs.append(("minced_rc", dr, 0))

        # Deduplicate by sequence
        seen = set()
        unique_drs = []
        for source, dr, n_rep in candidate_drs:
            if dr not in seen:
                seen.add(dr)
                unique_drs.append((source, dr, n_rep))

        log.info(f"  Found {len(unique_drs)} unique candidate DRs")

        protein_seq = _load_protein_seq(cid)
        scored = []
        for source, dr, n_rep in unique_drs:
            tag = f"T1_{source}_{len(dr)}nt"
            s = score_dr_with_protenix(cid, protein_seq, dr, info["subtype"], tag, skip_fold)
            s["source"] = f"tier1_{source}"
            s["n_repeats"] = n_rep
            scored.append(s)
            log.info(f"  DR {dr[:20]}... ({len(dr)}nt, {source}): ipTM={s['iptm']:.3f}")

        results[cid] = scored

    return results


# =====================================================================
# Tier 2: Salvage from original mining repeat_domains
# =====================================================================

def _cluster_overlapping_kmers(kmers: list) -> list:
    """Cluster overlapping sliding-window k-mers to reconstruct consensus repeats."""
    if not kmers:
        return []

    by_len = defaultdict(list)
    for k in kmers:
        k_clean = k.strip().upper()
        if k_clean and "N" not in k_clean:
            by_len[len(k_clean)].append(k_clean)

    consensus_candidates = []
    for klen, group in sorted(by_len.items()):
        if len(group) < 3:
            continue

        counts = Counter(group)
        if counts.most_common(1)[0][1] >= 2:
            top_kmer = counts.most_common(1)[0][0]
            consensus_candidates.append(top_kmer)
            continue

        # Try to reconstruct longer sequence from overlapping windows
        sorted_group = sorted(set(group))
        for seed in sorted_group[:5]:
            extended = seed
            for other in sorted_group:
                if other == seed:
                    continue
                overlap = min(len(seed), len(other)) - 1
                while overlap > 10:
                    if extended.endswith(other[:overlap]):
                        extended = extended + other[overlap:]
                        break
                    overlap -= 1
            if len(extended) > len(seed) and 23 <= len(extended) <= 50:
                consensus_candidates.append(extended)

    return list(set(consensus_candidates))


def tier2_salvage_mining_repeats(contigs_path: str, skip_fold: bool = False) -> dict:
    """Tier 2: Extract and validate DRs from original mining CSV repeat_domains."""
    log.info("=" * 60)
    log.info("TIER 2: Salvage from mining repeat_domains")
    log.info("=" * 60)

    contig_to_candidate = {info["contig"]: cid for cid, info in CANDIDATES.items()}
    csv_dir = PROJECT_ROOT / "data" / "mined_hits"
    contig_repeats = defaultdict(list)

    for csv_file in csv_dir.glob("*.csv"):
        try:
            with open(csv_file, encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    acc = row.get("sra_accession", "")
                    repeats_raw = row.get("repeat_domains", "")
                    if acc in contig_to_candidate and repeats_raw:
                        kmers = [k for k in repeats_raw.split("|") if k.strip()]
                        contig_repeats[acc].extend(kmers)
        except Exception:
            continue

    results = {}
    for contig_id, cid in contig_to_candidate.items():
        kmers = contig_repeats.get(contig_id, [])
        log.info(f"{cid}: {len(kmers)} raw k-mers from mining CSVs")

        candidate_drs = []

        # Method 1: existing salvage logic
        salvaged = _salvage_from_old_kmers(kmers)
        if salvaged:
            candidate_drs.append(("salvaged", salvaged))

        # Method 2: cluster overlapping k-mers for consensus
        reconstructed = _cluster_overlapping_kmers(kmers)
        for seq in reconstructed:
            rna = seq.replace("T", "U")
            if not is_trna_like(rna) and has_crispr_structure(rna):
                oriented = _pick_best_dr_orientation(rna)
                candidate_drs.append(("clustered", oriented))

        # Method 3: filter individual k-mers by length preference (28-36nt)
        seen = set()
        for k in kmers:
            k_rna = k.strip().upper().replace("T", "U")
            if k_rna in seen or len(k_rna) < 28 or len(k_rna) > 36:
                continue
            seen.add(k_rna)
            if not is_trna_like(k_rna) and has_crispr_structure(k_rna):
                oriented = _pick_best_dr_orientation(k_rna)
                candidate_drs.append(("filtered_kmer", oriented))
                if len(candidate_drs) >= 8:
                    break

        # Deduplicate
        final = []
        dr_seen = set()
        for source, dr in candidate_drs:
            if dr not in dr_seen:
                dr_seen.add(dr)
                final.append((source, dr))

        log.info(f"  {len(final)} validated DR candidates")

        subtype = CANDIDATES[cid]["subtype"]
        protein_seq = _load_protein_seq(cid)
        scored = []
        for source, dr in final[:6]:
            tag = f"T2_{source}_{len(dr)}nt"
            s = score_dr_with_protenix(cid, protein_seq, dr, subtype, tag, skip_fold)
            s["source"] = f"tier2_{source}"
            scored.append(s)
            log.info(f"  DR {dr[:20]}... ({len(dr)}nt, {source}): ipTM={s['iptm']:.3f}")

        results[cid] = scored

    return results


# =====================================================================
# Tier 3: Known Cas13 DR library screen
# =====================================================================

def tier3_known_dr_library(skip_fold: bool = False) -> dict:
    """Tier 3: Screen all known Cas13 DRs against each candidate."""
    log.info("=" * 60)
    log.info("TIER 3: Known Cas13 DR library screen")
    log.info(f"  Library size: {len(KNOWN_DR_LIBRARY)} DRs x {len(CANDIDATES)} candidates = {len(KNOWN_DR_LIBRARY) * len(CANDIDATES)} combinations")
    log.info("=" * 60)

    results = {}
    for cid, info in CANDIDATES.items():
        protein_seq = _load_protein_seq(cid)
        subtype = info["subtype"]
        scored = []

        for dr_name, dr_dna in KNOWN_DR_LIBRARY.items():
            dr_rna = dr_dna.upper().replace("T", "U")
            tag = f"T3_{dr_name}"
            s = score_dr_with_protenix(cid, protein_seq, dr_rna, subtype, tag, skip_fold)
            s["source"] = f"tier3_{dr_name}"
            scored.append(s)
            log.info(f"  {cid} x {dr_name} ({len(dr_rna)}nt): ipTM={s['iptm']:.3f}")

        scored.sort(key=lambda x: x["iptm"], reverse=True)
        results[cid] = scored

    return results


# =====================================================================
# Tier 4: Broader metagenomic CRISPR search
# =====================================================================

def tier4_metagenome_search(contigs_path: str, skip_fold: bool = False) -> dict:
    """Tier 4: Run CRISPR detection on ALL contigs, match by proximity and sample."""
    log.info("=" * 60)
    log.info("TIER 4: Broader metagenomic CRISPR search")
    log.info("=" * 60)

    minced_bin = shutil.which("minced")
    if not minced_bin:
        log.warning("MinCED not found, using Python detector only")

    all_arrays = []

    # Run minced on full contig set
    if minced_bin and os.path.exists(contigs_path):
        out_txt = contigs_path + ".minced.txt"
        out_gff = contigs_path + ".minced.gff"
        try:
            log.info(f"Running MinCED on full contig set: {contigs_path}")
            subprocess.run(
                [minced_bin, "-minNR", "3", contigs_path, out_txt, out_gff],
                capture_output=True, text=True, timeout=600,
            )

            if os.path.exists(out_txt):
                current_contig = ""
                current_repeat = ""
                with open(out_txt) as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("Sequence '"):
                            current_contig = line.split("'")[1]
                        elif line.startswith("CRISPR"):
                            current_repeat = ""
                        elif line and not line.startswith("POSITION") and not line.startswith("-") and not line.startswith("Repeats") and not line.startswith("Time"):
                            cols = line.split()
                            try:
                                pos = int(cols[0])
                                repeat_seq = cols[1]
                                if len(repeat_seq) >= 23 and set(repeat_seq.upper()) <= {"A", "C", "G", "T"}:
                                    rna = repeat_seq.upper().replace("T", "U")
                                    all_arrays.append({
                                        "contig": current_contig,
                                        "dr_rna": rna,
                                        "position": pos,
                                    })
                            except (ValueError, IndexError):
                                continue
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            log.warning(f"MinCED failed on full contig set: {e}")

    # Deduplicate DRs
    unique_drs = {}
    for arr in all_arrays:
        dr = arr["dr_rna"]
        if dr not in unique_drs:
            unique_drs[dr] = arr

    log.info(f"Found {len(unique_drs)} unique DRs across all contigs")

    # Filter to plausible Cas13 DRs (28-36nt, not tRNA, stem-loop)
    valid_drs = []
    for dr, info in unique_drs.items():
        if 23 <= len(dr) <= 50 and not is_trna_like(dr) and has_crispr_structure(dr):
            oriented = _pick_best_dr_orientation(dr)
            valid_drs.append((info["contig"], oriented))

    log.info(f"  {len(valid_drs)} pass structure filters")

    # Prefer DRs in Cas13 length range (28-36nt)
    preferred = [(c, d) for c, d in valid_drs if 28 <= len(d) <= 36]
    if preferred:
        valid_drs = preferred + [(c, d) for c, d in valid_drs if not (28 <= len(d) <= 36)]

    # Score top DRs per candidate (limit to 10 to manage compute)
    dr_set = list({d for _, d in valid_drs})[:10]
    log.info(f"  Scoring top {len(dr_set)} DRs against candidates")

    results = {}
    for cid, info in CANDIDATES.items():
        protein_seq = _load_protein_seq(cid)
        subtype = info["subtype"]
        scored = []

        for dr in dr_set:
            tag = f"T4_meta_{len(dr)}nt_{dr[:8]}"
            s = score_dr_with_protenix(cid, protein_seq, dr, subtype, tag, skip_fold)
            s["source"] = f"tier4_metagenome"
            scored.append(s)
            log.info(f"  {cid} x {dr[:20]}... ({len(dr)}nt): ipTM={s['iptm']:.3f}")

        scored.sort(key=lambda x: x["iptm"], reverse=True)
        results[cid] = scored

    return results


# =====================================================================
# Tier 5: Computational DR design
# =====================================================================

def _generate_stem_loop_candidates(stem_len: int = 6, loop_len: int = 4, flank_5: int = 4, flank_3: int = 4) -> list:
    """Generate RNA sequences with a single stem-loop structure."""
    bases = "ACGU"
    candidates = []

    stems = [
        ("GACCAC", "GUGGUC"),
        ("GAUUUA", "UAAAUC"),
        ("GUUUUG", "CAAAAC"),
        ("CCCCCA", "UGGGGG"),
        ("CAACCA", "UGGUUG"),
        ("GUUCAU", "AUGAAC"),
        ("GACCCC", "GGGGUC"),
        ("CCACCU", "AGGUGG"),
    ]
    loops = ["AAAA", "UUUU", "AAUG", "GAAA", "UGAA", "CCAA", "UAUC", "GCAG"]
    flanks_5 = ["GAUU", "AACC", "GUUU", "CAAC", "GCUA"]
    flanks_3 = ["AAAC", "UUAU", "GACC", "AACC", "CUAU"]

    for (stem5, stem3), loop, f5, f3 in product(stems, loops, flanks_5, flanks_3):
        dr = f5 + stem5 + loop + stem3 + f3
        if 26 <= len(dr) <= 38:
            candidates.append(dr)

    return candidates[:200]


def tier5_computational_design(best_per_candidate: dict, skip_fold: bool = False) -> dict:
    """Tier 5: Computational DR design starting from best known DR."""
    log.info("=" * 60)
    log.info("TIER 5: Computational DR design")
    log.info("=" * 60)

    results = {}

    for cid, info in CANDIDATES.items():
        protein_seq = _load_protein_seq(cid)
        subtype = info["subtype"]

        # Start from the best DR found so far
        best_so_far = best_per_candidate.get(cid, {})
        seed_dr = best_so_far.get("dr_rna", "")

        scored = []

        # Strategy 1: Mutate the best known DR
        if seed_dr and len(seed_dr) >= 20:
            log.info(f"{cid}: mutating seed DR ({len(seed_dr)}nt, ipTM={best_so_far.get('iptm', 0):.3f})")
            bases = "ACGU"
            mutants = set()
            for i in range(len(seed_dr)):
                for b in bases:
                    if b != seed_dr[i]:
                        mutant = seed_dr[:i] + b + seed_dr[i + 1:]
                        mutants.add(mutant)
            # Also try +-1nt truncations/extensions
            mutants.add(seed_dr[1:])
            mutants.add(seed_dr[:-1])
            for b in bases:
                mutants.add(b + seed_dr)
                mutants.add(seed_dr + b)

            # Filter mutants for structure
            valid_mutants = []
            for m in mutants:
                if 23 <= len(m) <= 50 and has_crispr_structure(m):
                    valid_mutants.append(m)

            # Score top 15 mutants (random sample if too many)
            import random
            if len(valid_mutants) > 15:
                valid_mutants = random.sample(valid_mutants, 15)

            for m in valid_mutants:
                tag = f"T5_mutant_{len(m)}nt"
                s = score_dr_with_protenix(cid, protein_seq, m, subtype, tag, skip_fold)
                s["source"] = "tier5_mutant"
                scored.append(s)

        # Strategy 2: De novo stem-loop generation
        log.info(f"{cid}: generating de novo stem-loop candidates")
        de_novo = _generate_stem_loop_candidates()
        valid_de_novo = [d for d in de_novo if has_crispr_structure(d)][:10]

        for dr in valid_de_novo:
            tag = f"T5_denovo_{len(dr)}nt"
            s = score_dr_with_protenix(cid, protein_seq, dr, subtype, tag, skip_fold)
            s["source"] = "tier5_denovo"
            scored.append(s)

        scored.sort(key=lambda x: x["iptm"], reverse=True)
        results[cid] = scored

        if scored and scored[0]["iptm"] > 0:
            log.info(f"  Best T5: ipTM={scored[0]['iptm']:.3f} ({scored[0]['dr_rna'][:20]}...)")

    return results


# =====================================================================
# Output: Report + metadata update + final base fold
# =====================================================================

def _write_report(all_results: list):
    """Write comprehensive CSV report of all tested DRs."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "candidate_id", "source", "tag", "dr_rna", "dr_len",
        "iptm", "ptm", "af2_ig", "json_path", "summary_path",
    ]
    with open(REPORT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in sorted(all_results, key=lambda x: (-x["iptm"], x["candidate_id"])):
            writer.writerow(row)
    log.info(f"Report written: {REPORT_CSV} ({len(all_results)} entries)")


def _update_metadata(best_per_candidate: dict):
    """Update variant_domain_metadata.json with best DRs."""
    if not METADATA_JSON.exists():
        log.warning(f"No metadata JSON at {METADATA_JSON}")
        return

    with open(METADATA_JSON) as f:
        metadata = json.load(f)

    updated = 0
    for cid, best in best_per_candidate.items():
        if cid in metadata and best.get("dr_rna"):
            metadata[cid]["crRNA_DR"] = best["dr_rna"]
            metadata[cid]["crRNA_repeat_used"] = best["dr_rna"]
            metadata[cid]["crrna_discovery_source"] = best.get("source", "unknown")
            metadata[cid]["crrna_discovery_iptm"] = best.get("iptm", 0.0)
            updated += 1

    with open(METADATA_JSON, "w") as f:
        json.dump(metadata, f, indent=2)
    log.info(f"Updated {updated} entries in {METADATA_JSON}")


def _update_db_crrna(best_per_candidate: dict):
    """Update SQLite DB with best crRNA repeats so JSON generation works."""
    if not DB_FILE.exists():
        return
    conn = sqlite3.connect(str(DB_FILE))
    cur = conn.cursor()
    for cid, best in best_per_candidate.items():
        if best.get("dr_rna"):
            dr_dna = best["dr_rna"].replace("U", "T")
            cur.execute("UPDATE variants SET crrna_repeat=? WHERE sequence_id=?", (dr_dna, cid))
    conn.commit()
    conn.close()
    log.info("Updated crRNA repeats in SQLite DB")


def _create_final_jsons(best_per_candidate: dict):
    """Create Protenix JSONs in jsons/ for the orchestrator pipeline."""
    JSONS_DIR.mkdir(parents=True, exist_ok=True)

    for cid, best in best_per_candidate.items():
        if not best.get("dr_rna"):
            continue

        protein_seq = _load_protein_seq(cid)
        subtype = CANDIDATES[cid]["subtype"]
        spacer = get_spacer_for_subtype(subtype)
        target_rna = get_target_for_spacer(spacer)
        crrna = assemble_crrna(best["dr_rna"], spacer, subtype)

        payload = [
            {
                "name": cid,
                "sequences": [
                    {"proteinChain": {"sequence": protein_seq, "count": 1}},
                    {"rnaSequence": {"sequence": crrna, "count": 1}},
                    {"rnaSequence": {"sequence": target_rna, "count": 1}},
                ],
            }
        ]
        out_path = JSONS_DIR / f"{cid}.json"
        with open(out_path, "w") as f:
            json.dump(payload, f, indent=2)
        log.info(f"Created pipeline JSON: {out_path}")


def _run_final_base_fold(best_per_candidate: dict):
    """Trigger final base-model Protenix fold for each candidate with best DR."""
    log.info("=" * 60)
    log.info("FINAL: Base-model structural validation with best DRs")
    log.info("=" * 60)

    for cid, best in best_per_candidate.items():
        if not best.get("dr_rna"):
            log.warning(f"{cid}: no DR found, skipping final fold")
            continue

        protein_seq = _load_protein_seq(cid)
        subtype = CANDIDATES[cid]["subtype"]
        spacer = get_spacer_for_subtype(subtype)
        crrna = assemble_crrna(best["dr_rna"], spacer, subtype)

        safe_tag = re.sub(r"[^\w\-.]", "_", f"{cid}_FINAL")
        json_path = STRUCTURE_CHECK_DIR / f"{safe_tag}.json"
        json_path.parent.mkdir(parents=True, exist_ok=True)

        payload = [
            {
                "name": safe_tag,
                "sequences": [
                    {"proteinChain": {"sequence": protein_seq, "count": 1}},
                    {"rnaSequence": {"sequence": crrna, "count": 1}},
                ],
            }
        ]
        with open(json_path, "w") as f:
            json.dump(payload, f, indent=2)

        log.info(f"Folding {cid} with best DR ({len(best['dr_rna'])}nt, source={best.get('source', '?')}, mini ipTM={best.get('iptm', 0):.3f})")
        try:
            out_dir = str(STRUCTURE_CHECK_DIR / f"{safe_tag}_results")
            struct, summary = run_protenix_inference(
                json_path=str(json_path),
                out_dir=out_dir,
                model_tier="base",
            )
            scores = extract_protenix_scores(summary)
            log.info(f"  FINAL {cid}: ipTM={scores['iptm']:.3f}  pTM={scores['ptm']:.3f}  AF2-IG={scores['af2_ig']:.3f}")
        except Exception as e:
            log.warning(f"  Final fold failed for {cid}: {e}")


# =====================================================================
# Main orchestrator
# =====================================================================

def main():
    parser = argparse.ArgumentParser(description="Iterative crRNA DR discovery for novel Cas13 candidates")
    parser.add_argument("--contigs", default="/tmp/all_contigs.fasta", help="Path to all contigs FASTA")
    parser.add_argument("--skip-fold", action="store_true", help="Skip Protenix folding (DR discovery only)")
    parser.add_argument("--tier", type=int, default=1, help="Start from this tier (1-5)")
    parser.add_argument("--threshold", type=float, default=IPTM_THRESHOLD, help=f"ipTM threshold for early stopping (default {IPTM_THRESHOLD})")
    parser.add_argument("--no-final-fold", action="store_true", help="Skip final base-model fold")
    args = parser.parse_args()

    threshold = args.threshold
    all_results = []
    best_per_candidate = {cid: {"iptm": 0.0, "dr_rna": ""} for cid in CANDIDATES}

    def _update_best(tier_results):
        """Merge tier results into global best and all_results."""
        for cid, scored in tier_results.items():
            all_results.extend(scored)
            for s in scored:
                if s["iptm"] > best_per_candidate[cid].get("iptm", 0):
                    best_per_candidate[cid] = s

    def _all_above_threshold():
        return all(best_per_candidate[cid]["iptm"] >= threshold for cid in CANDIDATES)

    def _report_status():
        log.info("-" * 50)
        log.info("Current best DRs:")
        for cid in CANDIDATES:
            b = best_per_candidate[cid]
            log.info(f"  {cid}: ipTM={b.get('iptm', 0):.3f} DR={b.get('dr_rna', 'none')[:25]}... ({b.get('source', '?')})")
        solved = sum(1 for b in best_per_candidate.values() if b.get("iptm", 0) >= threshold)
        log.info(f"  {solved}/{len(CANDIDATES)} above ipTM threshold ({threshold})")
        log.info("-" * 50)

    # Load contigs for tiers that need them
    contigs = {}
    if os.path.exists(args.contigs):
        log.info(f"Loading contigs from {args.contigs}")
        contigs = _load_contigs(args.contigs)
        log.info(f"Loaded {len(contigs)} contigs")
    else:
        log.warning(f"Contigs file not found: {args.contigs}")

    # --- Tier 1 ---
    if args.tier <= 1:
        results = tier1_proximity_detection(contigs, args.skip_fold)
        _update_best(results)
        _report_status()
        if _all_above_threshold():
            log.info("All candidates solved at Tier 1!")

    # --- Tier 2 ---
    if args.tier <= 2 and not _all_above_threshold():
        results = tier2_salvage_mining_repeats(args.contigs, args.skip_fold)
        _update_best(results)
        _report_status()
        if _all_above_threshold():
            log.info("All candidates solved at Tier 2!")

    # --- Tier 3 ---
    if args.tier <= 3 and not _all_above_threshold():
        results = tier3_known_dr_library(args.skip_fold)
        _update_best(results)
        _report_status()
        if _all_above_threshold():
            log.info("All candidates solved at Tier 3!")

    # --- Tier 4 ---
    if args.tier <= 4 and not _all_above_threshold():
        if contigs:
            results = tier4_metagenome_search(args.contigs, args.skip_fold)
            _update_best(results)
            _report_status()
            if _all_above_threshold():
                log.info("All candidates solved at Tier 4!")
        else:
            log.warning("Skipping Tier 4: no contigs available")

    # --- Tier 5 ---
    if args.tier <= 5 and not _all_above_threshold():
        results = tier5_computational_design(best_per_candidate, args.skip_fold)
        _update_best(results)
        _report_status()

    # --- Output ---
    log.info("=" * 60)
    log.info("DISCOVERY COMPLETE")
    log.info("=" * 60)

    _write_report(all_results)
    _update_metadata(best_per_candidate)
    _update_db_crrna(best_per_candidate)
    _create_final_jsons(best_per_candidate)

    for cid in CANDIDATES:
        b = best_per_candidate[cid]
        status = "SOLVED" if b.get("iptm", 0) >= threshold else "BEST_EFFORT"
        log.info(f"  [{status}] {cid}: ipTM={b.get('iptm', 0):.3f} DR={b.get('dr_rna', 'none')[:30]}... (source: {b.get('source', '?')})")

    # Final base-model fold
    if not args.no_final_fold and not args.skip_fold:
        _run_final_base_fold(best_per_candidate)


if __name__ == "__main__":
    main()
