<div align="center">

# CasCADE

### *The Dormant Blade Protocol*

**Cas** **C**ollateral **A**ctivation — **D**iscovery & **E**ngineering

*Weaponizing Cas13 collateral cleavage — an AI-driven structural pipeline for engineered suicide switches in targeted oncology.*

[![Rust](https://img.shields.io/badge/Rust-Powered-000000?logo=rust&logoColor=white)](https://github.com/twaite11/cattle-prod)
[![Cattle-Prod](https://img.shields.io/badge/Eval_Engine-Cattle--Prod-ff6b35)](https://github.com/twaite11/cattle-prod)
[![PXDesign](https://img.shields.io/badge/Generation-PXDesign-green)](https://github.com/bytedance/PXDesign)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue?logo=python&logoColor=white)](https://www.python.org)
[![Tests](https://img.shields.io/badge/tests-56%20passed-brightgreen?logo=pytest&logoColor=white)](#testing)
[![License](https://img.shields.io/badge/license-Research-lightgrey)](#license)

<br/>

[**Quick Start**](#-quick-start) · [**Architecture**](#-pipeline-architecture) · [**Deploy on VPS**](VPS_DEPLOY.md) · [**Interactive Diagram**](workflow_diagram.html) · [**RL Training Data**](RL_TRAINING_FORMAT.md)

<br/>

<p align="center"><sub>Metagenomic discovery → structural bootstrapping → RL evolution · Rust eval · SQLite dashboard</sub></p>

</div>

---

## 💡 The Core Insight

> Traditional gene therapies view Cas13's **collateral trans-cleavage** — the indiscriminate shredding of bystander RNA upon activation — as a catastrophic flaw. **CASCADE reverses this paradigm.**

We discover, validate, and engineer novel Cas13e-like proteins from metagenomic dark matter to function as **highly specific biological suicide switches**:

| State | HEPN Domains (His NE2-NE2) | Behavior |
|:------|:----------------------------|:---------|
| 🟢 **Healthy cell** (OFF) | **≥ 18 Å apart** — catalytically inert | Completely dormant. Zero leakiness. |
| 🔴 **Tumor cell** (ON) | **≤ 12 Å apart** — catalytically aligned | Massive collateral cleavage → programmed cell death |

The switch is triggered **only** when the engineered Cas13 detects and binds a **tumor-specific fusion RNA** (e.g., BCR-ABL, EWS-FLI1). The conformational change snaps the HEPN domains together, transforming a dormant ribonucleoprotein into a lethal RNA-shredding machine — strictly within the tumor.

---

## 🔬 Engineering Philosophy

<table>
<tr>
<td width="50%">

### 🔒 Nature Provides the Targeting
**Do Not Touch**

The wild-type crRNA binding pocket (REC lobe: NTD + Helical-1) was optimized by millions of years of evolution to secure its specific CRISPR repeat. During directed evolution, we **freeze the REC lobe entirely**. We use the native crRNA discovered from metagenomic contigs — no guide RNA redesign needed.

</td>
<td width="50%">

### ⚡ Physics Provides the Switch
**Engineer This**

The thermodynamic OFF→ON transition is controlled by **Inter-Domain Linkers (IDLs)** and the **Helical-2 domain**. We engineer these regions to:
- **Hyper-stabilize the dormant state** (high energetic barrier → zero leakiness)
- **Maximize trans-cleavage post-activation** (persistently locked in lethal state after tumor RNA forces the switch open)

</td>
</tr>
</table>

---

## 🏗️ Pipeline Architecture

CASCADE operates in three phases. Phase 0 discovers novel Cas13 enzymes from metagenomic dark matter. Phase 1 bootstraps validated baselines with native crRNAs and structural screening. Phase 2 runs an autonomous **reinforcement-learning (RL) driven evolution loop** that continuously designs, evaluates, and learns from structural predictions.

Structural evaluation is powered by **[Cattle-Prod](https://github.com/twaite11/cattle-prod)** — a Rust-native reimplementation of Protenix that compiles to a single binary with zero Python runtime overhead. The compound speedup on CPU-bound stages (parsing, tokenization, featurization, scoring) means more generations explored per GPU-hour.

```
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │  PHASE 0: Metagenomic Discovery                                             │
 │                                                                             │
 │  Raw Contigs (NCBI/SRA)                                                     │
 │       │                                                                     │
 │       ▼                                                                     │
 │  ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐              │
 │  │ DIAMOND  │───▶│   ORF    │───▶│ DIAMOND  │───▶│   HMM    │              │
 │  │ BLASTX   │    │ Extract  │    │ BLASTP   │    │ Classify │              │
 │  │ (6-frame)│    │ (±strand)│    │ (confirm)│    │ (subtype)│              │
 │  └──────────┘    └──────────┘    └──────────┘    └────┬─────┘              │
 │                                                       │                     │
 │       ┌───────────────────────────────────────────────┘                     │
 │       ▼                                                                     │
 │  ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐              │
 │  │  MinCED  │───▶│    DR    │───▶│ Protein  │───▶│ Dedup +  │              │
 │  │  CRISPR  │    │ Assign + │    │ Validate │    │ Rank by  │              │
 │  │  Arrays  │    │ 5-Tier   │    │ (QC +    │    │ Confid-  │              │
 │  │ Detection│    │ crRNA    │    │  HEPN    │    │ ence     │──▶ Phase 1   │
 │  │          │    │ Discovery│    │  Archit.)│    │ (H/M/L)  │              │
 │  └──────────┘    └──────────┘    └──────────┘    └──────────┘              │
 │                                                                             │
 │  Confidence tiers:                                                          │
 │    HIGH   = Cas13 BLAST hit + adjacent CRISPR array + dual HEPN motifs     │
 │    MEDIUM = Cas13 BLAST hit + dual HEPN (no adjacent array)                │
 │    LOW    = Dual HEPN only (exploratory)                                    │
 └──────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │  PHASE 1: Bootstrap Baselines                                               │
 │                                                                             │
 │  ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐              │
 │  │ Parse &  │───▶│ HMM      │───▶│ Fix crRNA│───▶│ Structur-│              │
 │  │ Annotate │    │ Classify │    │ Assign-  │    │ al Screen│              │
 │  │ (SQLite +│    │ (subtype │    │ ments    │    │ (mini    │              │
 │  │  HEPN    │    │  + HEPN  │    │ (offline/│    │  fold)   │              │
 │  │  anchor) │    │  PF05168)│    │  fetch)  │    │          │              │
 │  └──────────┘    └──────────┘    └──────────┘    └────┬─────┘              │
 │                                                       │                     │
 │       ┌───────────────────────────────────────────────┘                     │
 │       ▼                                                                     │
 │  ┌──────────┐    ┌──────────┐                                               │
 │  │ MSA Re-  │───▶│ Validate │                                               │
 │  │ run Top  │    │ CRISPR   │──▶ Validated baselines                        │
 │  │ N (base  │    │ Repeats  │    with native crRNAs                         │
 │  │  model)  │    │ (RNAfold)│    + Phase 1 structures                       │
 │  └──────────┘    └──────────┘                                               │
 └──────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │  PHASE 2: Active Learning Evolution Loop                                    │
 │                                                                             │
 │  For each lineage (parallel workers with GPU lock):                         │
 │                                                                             │
 │  ┌────────────────────────────────────────────────────────────┐             │
 │  │ Gen 0: Baseline Reference                                  │             │
 │  │   Evaluate unmodified enzyme → set RL fitness baseline     │             │
 │  └────────────────────┬───────────────────────────────────────┘             │
 │                       ▼                                                     │
 │  ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐             │
 │  │ PXDesign │───▶│  HEPN    │───▶│  RL Bias │───▶│ Protenix │             │
 │  │ Backbone │    │ Stitch + │    │  Apply   │    │ Mini OFF │             │
 │  │(Gen 1) or│    │ Resolve  │    │ (linker  │    │ + ON     │             │
 │  │ MPNN Re- │    │ X → G/WT │    │  regions │    │ eval     │             │
 │  │fine(Gen2+│    │          │    │  only)   │    │          │             │
 │  └──────────┘    └──────────┘    └──────────┘    └────┬─────┘             │
 │                                                       │                    │
 │       ┌───────────────────────────────────────────────┘                    │
 │       ▼                                                                    │
 │  ┌──────────┐         ┌──────────┐    ┌──────────┐    ┌──────────┐        │
 │  │ HEPN 3D  │──pass──▶│ Protenix │───▶│ Off-tgt  │───▶│ Fitness  │        │
 │  │ Distance │         │ Base     │    │ Specific.│    │ Score    │        │
 │  │ Filter   │         │ Ternary  │    │ (1/2/3mm │    │ (compos- │        │
 │  │ OFF≥18Å  │         │ Complex  │    │ penalty) │    │  ite)    │        │
 │  │ ON ≤12Å  │         │ Eval     │    │          │    │          │        │
 │  └──────────┘         └──────────┘    └──────────┘    └────┬─────┘        │
 │       │ fail                                               │               │
 │       ▼                                                    ▼               │
 │  Score with                                          ┌──────────┐          │
 │  mini metrics                                        │Evolution │          │
 │  (no base eval)────────────────────────────────────▶│  Gym RL  │          │
 │                                                      │ (relative│          │
 │  ┌──────────┐    ┌──────────┐    ┌──────────┐       │  fitness │          │
 │  │Population│◀───│ Global   │◀───│ RL Bias  │◀──────│  weights)│          │
 │  │ Update   │    │ Best     │    │ Matrix   │       └──────────┘          │
 │  │ (top-K   │    │ Tracking │    │ Export   │                              │
 │  │ tourney) │    │ (never   │    │ (PSSM →  │                              │
 │  │          │    │ regress) │    │  MPNN)   │                              │
 │  └────┬─────┘    └──────────┘    └──────────┘                              │
 │       │                                                                    │
 │       ▼                                                                    │
 │  ┌──────────────────────────────────────────────┐                          │
 │  │ Stagnation? ──yes──▶ Abandon lineage         │                          │
 │  │ Elite found? ──yes──▶ Save + next lineage     │                          │
 │  │ Otherwise ──────────▶ Next generation (loop)  │                          │
 │  └──────────────────────────────────────────────┘                          │
 └──────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
                          ┌──────────────────────┐
                          │  outputs/             │
                          │  optimized_switches/  │
                          │   *_optimal.fasta     │
                          │   *_ternary.cif       │
                          │   *_crRNA.fasta       │
                          │  rl_gym_data/          │
                          │   rl_training.jsonl    │
                          │  switch_report.csv     │
                          └──────────────────────┘
```

> 📊 **[Open the interactive workflow diagram →](workflow_diagram.html)** for the full node-by-node breakdown with color-coded phases.

### Phase 0 — Metagenomic Discovery

| Step | Script | What It Does |
|:-----|:-------|:-------------|
| **0a** | `mining_v2.py` | Full discovery pipeline: DIAMOND BLASTX on raw contigs → 6-frame ORF extraction (both strands) → HMM classification → MinCED CRISPR array detection → DR assignment → confidence-ranked hits. |
| **0b** | `blast_candidates.py` | DIAMOND BLASTP confirmation of extracted ORFs against comprehensive Cas13 reference DB (8 characterized proteins across Cas13a/b/d/bt subtypes). |
| **0c** | `classify_cas13_hmm.py` | HMM domain classification using Pfam HEPN (PF05168) + Cas13 family-specific diagnostic motifs. Assigns subtype (a/b/d/X/Y) and validates dual-HEPN architecture. |
| **0d** | `validate_proteins.py` | Protein QC: length range, composition, low-complexity filter, HEPN motif spacing (150-600 residues), alpha-helical propensity, charge profile. |
| **0e** | `discover_crrna.py` | 5-tier iterative crRNA DR discovery: (1) proximity CRISPR detection, (2) mining DR salvage, (3) known DR library screen via Protenix ipTM, (4) broader metagenome search, (5) computational stem-loop design + mutagenesis. |
| **0f** | `fix_crrna_assignments.py` | Replace naive k-mer "repeats" (often tRNAs) with properly detected CRISPR direct repeats. Fetches source contigs, detects arrays, filters tRNA-like sequences, re-assigns DRs by proximity. |

### Phase 1 — Bootstrap Baselines

| Step | Script | What It Does |
|:-----|:-------|:-------------|
| **1a** | `01_parse_and_annotate.py` | Parse FASTA + CSV into SQLite DB. Anchor HEPN1/HEPN2 domains via `R.{4,6}H` motif with `_select_hepn_pair()` (positional + separation constraints, not naive first/last). Generate prediction-compatible JSONs with protein + crRNA + target RNA. |
| **1b** | `classify_cas13_hmm.py` | *(If not run in Phase 0)* HMM subtype classification and HEPN domain validation. |
| **1c** | `fix_crrna_assignments.py --offline` | Validate/repair crRNA DR assignments using improved filters (tRNA exclusion, array structure). Updates `variant_domain_metadata.json`. |
| **1d** | `02_run_screening.sh` | GPU-accelerated Cattle-Prod/Protenix mini structural screen. Filters hits that can't form bilobed structures or bind crRNA. Generates Phase 1 structures in `outputs/phase1_screening/`. |
| **1e** | `02b_rerun_top_with_msa.sh` | *(Optional)* Re-run top N baselines with MSA-enhanced base model for higher-quality seed structures and more reliable ipTM/pTM scores. |
| **1f** | `validate_crispr_repeats.py` | *(Optional)* Validate CRISPR repeats via RNAfold stem-loop structure. Outputs `validated_baseline_ids.txt` to restrict evolution to structurally validated DRs. |

### Phase 2 — Active Learning Evolution Loop

The evolution orchestrator (`evolution_orchestrator.py`) runs an autonomous loop with parallel workers, tournament selection, and closed-loop RL:

```
Gen 0 — Baseline Reference:
  1. EVALUATE  → Protenix mini predicts unmodified enzyme OFF + ON structures
  2. MEASURE   → 3D HEPN NE2-NE2 distance + ipTM + AF2-IG scores
  3. BASELINE  → Set RL fitness reference (EvolutionGym centers rewards here)

For each generation (1 → MAX_GENERATIONS):
  1. SELECT    → Tournament selection picks parent from population (top-K pool)
  2. GENERATE  → Gen 1: PXDesign diffusion → novel backbone geometry
                 Gen 2+: MPNN iterative refinement on cached backbone
                          (exploration schedule unfreezes 10-50% of linker positions)
  3. STITCH    → Wild-type HEPN1/HEPN2 catalytic domains grafted into designed linkers
  4. RESOLVE   → Unknown residues (X) replaced from baseline or Glycine fallback
  5. BIAS      → RL bias matrix applied as PSSM to linker regions (closed-loop RL)
  6. EVALUATE  → Protenix mini predicts OFF-state and ON-state structures
  7. MEASURE   → 3D HEPN His NE2 distance (Å) + ipTM + AF2-IG confidence scores
  8. FILTER    → OFF ≥ 18Å (dormant) AND ON ≤ 12Å (active) → passes to base eval
  9. PROMOTE   → Passing variants: Protenix base ternary complex (protein+crRNA+target)
  10. TEST     → Off-target specificity: 1/2/3-mismatch guides (progressive penalty)
  11. SCORE    → Composite fitness: HEPN_shift + ipTM×50 + AF2-IG×20 − penalties
  12. LEARN    → EvolutionGym: relative fitness normalization → mutation weight update
  13. EXPORT   → Bias matrix → mpnn_bias_gen_X.json (position→AA PSSM for MPNN)
  14. ADVANCE  → Population merge (top-K survive). Global best → next baseline.
  15. CHECK    → Elite found (ipTM≥0.85, AF2-IG≥0.80, ON≤12Å)? → save + next lineage
                 Stagnated (4 gens no improvement)? → abandon lineage
```

<details>
<summary><b>🧮 Fitness Function Details</b></summary>

```
fitness = (off_dist − on_dist) − (MIN_OFF − MAX_ON)     # HEPN conformational shift
        + (ipTM − 0.7) × 50                              # Structure prediction confidence
        + (af2_ig − 0.5) × 20                            # Interface confidence
        × 2.0 if full ternary evaluation                  # Bonus for complete complex
        − specificity_penalties                            # Progressive: 3mm > 2mm > 1mm
        − 5.0 if fallback variant                          # Penalty for stitch failures
```

**Thresholds (NE2-NE2 measurement):**
- OFF distance ≥ 18.0 Å (catalytically inert when dormant)
- ON distance ≤ 12.0 Å (catalytically aligned when activated)
- Elite: ipTM ≥ 0.85, AF2-IG ≥ 0.80, ON ≤ 12.0 Å

</details>

<details>
<summary><b>🔄 RL Bias Loop Details</b></summary>

The `EvolutionGym` tracks every mutation across all generations:

1. **Beneficial mutations** (from high-fitness variants) receive positive weight
2. **Harmful mutations** (from low-fitness variants) receive negative weight
3. Weights are exported as `mpnn_bias_gen_X.json` — a position→amino-acid bias matrix
4. The bias is applied **post-stitching** to linker regions only (HEPN catalytic sites are never modified)
5. Only substitution mutations are tracked (insertions/deletions are excluded from the bias matrix)
6. Bias threshold: weight > 0.5 required before applying a substitution

This creates a **closed-loop reinforcement signal**: good designs → stronger bias toward their mutations → better designs.

</details>

<details>
<summary><b>🏆 Global Best Tracking</b></summary>

The pipeline tracks the **all-time highest-fitness protein** across all generations and lineages:

- After each generation, the best variant is compared against the global best
- If the new variant is better → it becomes the new global best and next baseline
- If not → the global best is reused as the next baseline (never regresses)
- Elite variants (ipTM ≥ 0.85, AF2-IG ≥ 0.80, ON ≤ 12 Å) are saved to `optimized_switches/`
- On lineage exhaustion, the pipeline falls back to the global best before trying new lineages

</details>

---

## 🚀 Quick Start

### Prerequisites

| Requirement | Minimum |
|:------------|:--------|
| **GPU** | A100 40GB (80GB recommended) |
| **OS** | Ubuntu 22.04+ with CUDA 12.x |
| **Rust** | 1.78+ (for Cattle-Prod) |
| **Python** | 3.11+ |
| **Conda** | Miniconda / Anaconda |
| **Storage** | 100–200 GB free |

### Setup (One-Time)

```bash
git clone https://github.com/twaite11/CASCADE-Cas-Collateral-Activation-Discovery-Engineering.git
cd CASCADE-Cas-Collateral-Activation-Discovery-Engineering

# 1. Build Cattle-Prod (Rust eval engine)
git clone https://github.com/twaite11/cattle-prod.git ../cattle-prod
cd ../cattle-prod/cattle-prod && cargo build --release -p cattle-prod-cli && cd -
export PATH="$(realpath ../cattle-prod/cattle-prod/target/release):$PATH"

# 2. Automated dual-environment setup (detects Cattle-Prod + installs PXDesign)
chmod +x scripts/setup_dual_env.sh
./scripts/setup_dual_env.sh
```

The setup script auto-detects Cattle-Prod and configures two conda environments:

| Component | Language | Purpose |
|:----------|:---------|:--------|
| **Cattle-Prod** | Rust | Structure prediction engine (single binary, ~5ms startup) |
| `cascade` env | Python | Orchestration, fitness scoring, evolution loop |
| `pxdesign` env | Python | PXDesign variant generation (called cross-env) |

> **Why Cattle-Prod over Protenix?** Cattle-Prod is a Rust reimplementation of Protenix that eliminates Python interpreter overhead, parallelizes CPU-bound featurization via Rayon, and deploys as a single 15MB binary. CASCADE auto-detects it on PATH and falls back to Protenix if unavailable.

### Activate (Every Session)

```bash
source scripts/cascade_env.sh
```

### Run

```bash
cd scripts

# ── Phase 0: Discover novel Cas13 from metagenomic contigs ──
# Download contigs into data/new_contigs/, then:
python mining_v2.py --contigs ../data/new_contigs/*.fasta --output-dir ../outputs/mining_v2
python blast_candidates.py                        # BLASTP confirmation
python classify_cas13_hmm.py                      # HMM subtype assignment
python validate_proteins.py                       # Protein QC
# Place validated hits in data/mined_hits/ (FASTA + metadata CSV)

# ── Phase 1: Bootstrap validated baselines ──
python 01_parse_and_annotate.py                   # SQLite + HEPN anchoring + JSONs
python classify_cas13_hmm.py                      # Subtype classification (if not done above)
python fix_crrna_assignments.py --offline          # Repair crRNA DR assignments
./02_run_screening.sh                             # Structural screen (mini model)
./02b_rerun_top_with_msa.sh 4                     # MSA re-run top 4 (optional)

# ── Phase 2: Active learning evolution ──
CASCADE_WORKERS=1 python evolution_orchestrator.py
```

Or run Phases 1+2 at once:

```bash
./scripts/run_pipeline.sh 2>&1 | tee "logs/cascade_$(date +%Y%m%d_%H%M%S).log"
```

> 📖 **Full deployment guide:** [VPS_DEPLOY.md](VPS_DEPLOY.md)

---

## Live Dashboard (SQLite-First Architecture)

CASCADE now includes a thin dashboard backend designed for VPS operation:

- Pipeline writes continue to `SQLite` (`metadata/cas13_variants.db`)
- `FastAPI` service reads SQLite in `WAL` mode and serves aggregated endpoints
- Web frontend auto-refreshes every N seconds for near-real-time visibility
- Storage is adapter-based so `Postgres` can be added later without rewriting API handlers
- Service is read-only (no mutation of pipeline outputs)

Start the dashboard:

```bash
pip install -r requirements.txt
uvicorn dashboard_backend.main:app --host 0.0.0.0 --port 8000
```

Then open:

- API health: `http://<server-ip>:8000/health`
- Dashboard UI: `http://<server-ip>:8000/dashboard`
- Runbook: `dashboard_backend/README.md`

### Containerized controller + Vast.ai parallel runs

For one-click bring-up of the full controller stack (dashboard + Vast.ai
orchestration + React SPA) and to launch each parallel evolution run on a
dedicated A100 VPS, use the bundled Docker image:

```bash
cp .env.example .env            # edit VAST_SSH_KEY_PATH, etc.
docker compose up -d --build    # http://localhost:8000/dashboard
```

From the dashboard you can:

1. **Baselines tab** — pick enzyme IDs (and their matching crRNA) to evolve.
2. **Launch runs** — the controller calls `vastai create instance` per run,
   streams orchestrator logs via WebSocket (xterm.js), and rsyncs results back
   on completion. No Redis, no Celery — SQLite + asyncio only.
3. **Optimized Switches sidebar** — live leaderboard with per-variant 3Dmol
   toggle (ON / OFF / off-target), HEPN1/HEPN2 coloring, and filters.
4. **Variant detail drawer** — full metadata, side-by-side 3D compare,
   crRNA spacer viz, artifact downloads. Press `?` for keyboard shortcuts.

> 📖 **Full container deployment guide:** [CONTAINER_DEPLOY.md](CONTAINER_DEPLOY.md)

### Concise Variant Naming

Generated variant IDs are now concise and generation-stable:

- `Lxxxxxx_g01_v00`
- `Lxxxxxx_g01_v00_fb` (fallback sequence)

This replaces compounding names like `..._variant_0_variant_1_ON` and keeps per-generation tracking readable.

---

## 📂 Project Structure

```
CASCADE/
│
├── 📄 README.md                        # You are here
├── 📄 VPS_DEPLOY.md                    # GPU VPS deployment guide
├── 📄 RL_TRAINING_FORMAT.md            # Post-training data format (DRAKES / ProteinMPNN)
├── 📄 workflow_diagram.html            # Interactive Mermaid workflow diagram
├── 📄 requirements.txt                 # Python dependencies
├── 📁 dashboard_backend/               # FastAPI aggregation service (SQLite WAL)
├── 📁 dashboard_frontend/              # Auto-refresh web UI
│
├── 📁 data/
│   ├── 📁 mined_hits/                  # Input: validated Cas13 FASTAs + metadata CSVs
│   │   ├── deep_hits_*.fasta           #   Protein ORFs from mining pipeline
│   │   └── deep_hits_*_metadata.csv    #   Accessions, DRs, confidence tiers
│   ├── 📁 hmm/                         # HMM profiles for domain classification
│   │   └── HEPN.hmm                    #   Pfam PF05168 HEPN domain profile
│   ├── 📁 new_contigs/                 # Raw metagenomic contigs (NCBI/SRA downloads)
│   └── cas13_reference_db.fasta        # Comprehensive Cas13 reference DB (8 proteins, 4 subtypes)
│
├── 📁 rust/                            # Optional Rust accelerators (cargo build --release)
│   ├── cascade_ingest/                 #   FASTA/CSV → SQLite + HEPN scanning + JSON gen
│   ├── cascade_sequtils/               #   Mutation extraction, motif finding, k-mer validation
│   └── cascade_structscore/            #   CIF/PDB parsing, CA-CA distance, score extraction
│
├── 📁 scripts/
│   ├── 🔧 setup_dual_env.sh           # Creates cascade + pxdesign conda envs
│   ├── 🔧 setup_vps.sh                # Simple single-env venv setup
│   ├── 🔧 cascade_env.sh              # Source this to activate env + PXDESIGN_CMD
│   ├── 🔧 run_pipeline.sh             # Full pipeline runner with logging
│   │
│   │── ── Phase 0: Metagenomic Discovery ──
│   ├── 📜 mining_v2.py                # Full discovery: BLASTX → ORF extract → HMM → CRISPR → rank
│   ├── 📜 blast_candidates.py          # DIAMOND BLASTP/BLASTX against Cas13 reference DB
│   ├── 📜 classify_cas13_hmm.py        # HMM subtype classification (HEPN PF05168 + family motifs)
│   ├── 📜 validate_proteins.py         # Protein QC: length, composition, HEPN architecture
│   ├── 📜 discover_crrna.py            # 5-tier crRNA DR discovery (proximity → library → design)
│   ├── 📜 fix_crrna_assignments.py     # Replace naive k-mers with real CRISPR DRs
│   │
│   │── ── Phase 1: Bootstrap Baselines ──
│   ├── 📜 01_parse_and_annotate.py     # Ingest → SQLite DB + HEPN anchoring + JSON generation
│   ├── 📜 02_run_screening.sh          # Cattle-Prod/Protenix mini structural screen
│   ├── 📜 02b_rerun_top_with_msa.sh   # MSA-enhanced base model re-run for top N
│   ├── 📜 validate_crispr_repeats.py   # CRISPR repeat validation via RNAfold
│   │
│   │── ── Phase 2: Evolution Loop ──
│   ├── 📜 03_pxdesign_wrapper.py       # PXDesign backbone + MPNN refinement + HEPN stitching
│   ├── 📜 evolution_orchestrator.py    # Active learning master controller (EvolutionGym RL)
│   ├── 📜 monitor_variants.py          # Live monitoring of evolution progress
│   │
│   └── 📁 utils/
│       ├── 📜 protenix_eval.py         # ON/OFF JSON payloads + Cattle-Prod/Protenix inference
│       ├── 📜 pdb_kinematics.py        # 3D HEPN NE2 distance + ipTM/pTM/AF2-IG extraction
│       └── 📜 hepn_structural_stitch.py # Graft WT HEPN catalytic domains into designed linkers
│
├── 📁 tests/                           # 56 unit tests (no GPU required)
│   ├── conftest.py                     # Shared fixtures
│   ├── test_01_parse_and_annotate.py
│   ├── test_evolution_orchestrator.py
│   ├── test_pxdesign_wrapper.py
│   ├── test_hepn_stitch.py
│   ├── test_pdb_kinematics.py
│   ├── test_protenix_eval.py
│   └── test_run_protenix_mocked.py
│
├── 📁 metadata/                        # Generated at runtime
│   ├── cas13_variants.db               #   SQLite database (all candidates + annotations)
│   └── variant_domain_metadata.json    #   HEPN domain boundaries + crRNA + subtype per variant
│
├── 📁 jsons/                           # Protenix/Cattle-Prod input payloads (protein+crRNA+target)
│
└── 📁 outputs/
    ├── mining_v2/                      # Phase 0: mining results, BLAST outputs, CRISPR arrays
    ├── phase1_screening/               # Phase 1: baseline CIF/PDB structures (mini model)
    ├── backbone_cache/                 # Cached PXDesign backbones for MPNN refinement (Gen 2+)
    ├── generation_queue/               # PXDesign/MPNN variant FASTAs per worker per generation
    ├── fast_eval/                      # Mini-model OFF/ON structural screening
    ├── high_fidelity_scoring/          # Base-model ternary complex predictions
    ├── rl_gym_data/                    # RL bias matrices + training dataset
    │   ├── mpnn_bias_gen_X.json        #   Per-generation PSSM bias matrices
    │   └── rl_training_dataset.jsonl   #   For DRAKES/ProteinMPNN post-training
    ├── optimized_switches/             # Elite outputs
    │   ├── *_optimal.fasta             #   Best engineered protein sequences
    │   ├── *_ternary_complex.cif       #   Predicted 3D ternary structures
    │   └── *_crRNA.fasta               #   Native crRNA sequences
    ├── switch_report.csv               # Full variant-level scoring report (all generations)
    ├── crrna_discovery_report.csv      # Phase 0: crRNA discovery results per tier
    └── cas13_classification_report.csv # Phase 0: HMM classification results
```

---

## ⚡ Rust-Powered Stack

CASCADE is built on a Rust-first philosophy. The most critical component — structure prediction — runs as a compiled Rust binary, and three additional Rust tools accelerate the CPU-bound orchestration. Python handles glue logic and PXDesign integration.

### Cattle-Prod (Eval Engine)

**[Cattle-Prod](https://github.com/twaite11/cattle-prod)** is the structure prediction engine. It's a full Rust reimplementation of Protenix (AlphaFold3-class) with:

| | Python (Protenix) | Rust (Cattle-Prod) |
|---|---|---|
| **Startup** | ~4s | ~5ms |
| **CPU featurization** | GIL-bound | Rayon parallel |
| **Deployment** | conda + 2GB deps | Single 15MB binary |
| **Memory safety** | Runtime errors | Compile-time guarantees |

CASCADE auto-detects `cattle-prod` on PATH. Override with `EVAL_CMD`:

```bash
export EVAL_CMD=/path/to/cattle-prod   # or: export EVAL_CMD=protenix (fallback)
```

### Pipeline Accelerators (Optional)

Three additional Rust CLI tools accelerate CPU-bound orchestration steps. When present on `PATH` (or in `rust/target/release/`), Python scripts automatically delegate to them. **Zero breaking changes.**

```bash
cd rust && cargo build --release
```

| Binary | Accelerates | Speedup |
|:-------|:------------|:--------|
| `cascade_ingest` | Phase 1a: FASTA/CSV parsing, HEPN regex scanning, SQLite bulk insert, JSON generation | 10-50x |
| `cascade_sequtils` | Evolution loop: mutation extraction, histidine motif finding, CRISPR repeat k-mer validation | 5-20x |
| `cascade_structscore` | Evolution loop: CIF/PDB structure parsing, CA-CA distance calculation, score extraction | 5-15x |

> `setup_dual_env.sh` automatically builds Cattle-Prod and these accelerators if `cargo` is installed.

---

## 🧪 Testing

All tests run on CPU — no GPU, Protenix, or PXDesign installation required.

```bash
pip install pytest
pytest tests/ -v
```

<details>
<summary><b>Test Coverage (56 tests)</b></summary>

| Module | Tests | What's Covered |
|:-------|:------|:---------------|
| `01_parse_and_annotate` | 6 | DB schema, FASTA/CSV loading, HEPN motif detection, JSON generation |
| `evolution_orchestrator` | 13 | Fitness computation, EvolutionGym EMA, bias clipping, mutation extraction, HEPN catalytic indices, metadata override |
| `03_pxdesign_wrapper` | 10 | Freeze config, HEPN boundary respect, X-residue resolution, RL bias application |
| `hepn_structural_stitch` | 2 | Domain grafting, short-sequence handling |
| `pdb_kinematics` | 7 | HEPN distance calculation, Protenix score extraction, AF2-IG fallback |
| `protenix_eval` | 10 | Mismatch generation, OFF/ON JSON payloads, off-target JSON, crRNA lookup |
| `run_protenix (mocked)` | 2 | CLI command construction for mini + base tiers, `pred`/`predict` fallback |

</details>

---

## 📊 Outputs

After the evolution loop completes, elite switches are saved to `outputs/optimized_switches/`:

| File | Description |
|:-----|:------------|
| `*_optimal.fasta` | Engineered Cas13 protein sequence — ready for synthesis |
| `*_ternary_complex.cif` | Predicted 3D structure (protein + crRNA + tumor RNA) |
| `*_crRNA.fasta` | Native crRNA sequence for this Cas13 variant |

The RL training dataset (`outputs/rl_gym_data/rl_training_dataset.jsonl`) can be used for **post-training** ProteinMPNN or discrete diffusion models like [DRAKES](https://github.com/ChenyuWang-Monica/DRAKES). See [RL_TRAINING_FORMAT.md](RL_TRAINING_FORMAT.md) for the schema.

---

## ⚙️ Key Configuration

These constants in `evolution_orchestrator.py` control the evolution:

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `MAX_GENERATIONS` | 12 | Generations per lineage (env: `CASCADE_MAX_GENERATIONS`) |
| `VARIANTS_PER_GEN` | 5 | Designs per generation (env: `CASCADE_VARIANTS_PER_GEN`) |
| `NUM_WORKERS` | 3 | Parallel lineage workers (env: `CASCADE_WORKERS`) |
| `POPULATION_SIZE` | 3 | Top-K survivors per generation |
| `TOURNAMENT_SIZE` | 2 | Parent selection pressure |
| `STAGNATION_LIMIT` | 4 | Gens without improvement before abandoning lineage |
| `MIN_OFF_DISTANCE` | 18.0 Å | Minimum HEPN NE2-NE2 distance for dormant state |
| `MAX_ON_DISTANCE` | 12.0 Å | Maximum HEPN NE2-NE2 distance for active state |
| `MIN_IPTM_SCORE` | 0.85 | Elite threshold for interface prediction confidence |
| `MIN_AF2_IG_SCORE` | 0.80 | Elite threshold for interface quality |
| `FALLBACK_FITNESS_PENALTY` | 5.0 | Penalty for stitch-failure fallback variants |
| `SPECIFICITY_PENALTY_BASE` | 0.3 | Off-target penalty (scaled by mismatch count) |

Environment variables:

| Variable | Example | Description |
|:---------|:--------|:------------|
| `CASCADE_WORKERS` | `1` | Number of parallel worker processes (use 1 on single GPU) |
| `CASCADE_MAX_GENERATIONS` | `12` | Override max generations per lineage |
| `CASCADE_VARIANTS_PER_GEN` | `3` | Override designs per generation |
| `EVAL_CMD` | `cattle-prod` or `/path/to/cattle-prod` | Structure prediction engine (auto-detected if on PATH) |
| `PXDESIGN_CMD` | `/path/to/envs/pxdesign/bin/pxdesign` | Direct path to PXDesign binary in its conda env |
| `PXDESIGN_SUBCOMMAND` | `pipeline` or `infer` | PXDesign CLI subcommand (auto-detected) |
| `PROTEINMPNN_DIR` | `/workspace/ProteinMPNN` | Path to ProteinMPNN repo (auto-detected) |
| `CATTLE_PROD_MINI_CKPT` | `/workspace/models/cattle/mini` | Cattle-prod mini checkpoint directory |
| `CATTLE_PROD_BASE_CKPT` | `/workspace/models/cattle/base` | Cattle-prod base checkpoint directory |
| `CATTLE_PROD_STRICT` | `1` | `1`: fail fast on errors. `0`: fallback to Protenix. |
| `CUDA_VERSION` | `12.1` | CUDA version for dual-env setup |

---

## Cattle-Prod Stability Runbook

Use this when rolling out new checkpoints:

```bash
# 1) Convert and diagnose
python /workspace/cattle-prod/cattle-prod/scripts/convert_weights.py model.pt \
  -o /workspace/models/cattle/base/model.safetensors \
  --mapping_version v1 \
  --diagnostic_report /workspace/models/cattle/base/diag.json

# 2) Verify converter output key coverage
python /workspace/cattle-prod/cattle-prod/scripts/convert_weights.py \
  --verify_safetensors /workspace/models/cattle/base/model.safetensors \
  --mapping_version v1

# 3) Verify runtime loadability from Rust side
cattle-prod verify-checkpoint \
  --checkpoint /workspace/models/cattle/base \
  -n cattle_prod_base_default_v1.0.0
```

If (3) fails:
- In production: set `CATTLE_PROD_STRICT=0` for explicit fallback while you remap the checkpoint.
- For hard-fail enforcement: keep `CATTLE_PROD_STRICT=1`.

Release gate before enabling strict cattle-prod mode:
- converter diagnostics show no missing required keys
- `verify-checkpoint` passes for mini/base checkpoint dirs
- one CASCADE generation produces non-default confidence metrics (not all zeros / constant fallback-style fitness)

---

## 🧬 Expected Input Data

### Option A: Start from raw contigs (recommended — full Phase 0)

Place assembled metagenomic contigs in `data/new_contigs/`:

```
data/new_contigs/
  SRR12345678_contigs.fasta      # Assembled contigs from SRA
  GCA_012345678_genomic.fna      # NCBI genome assemblies
  custom_metagenome.fasta        # Any assembled DNA sequences
```

Run `mining_v2.py` to discover Cas13 candidates automatically. The pipeline extracts ORFs, runs DIAMOND BLAST, detects CRISPR arrays, and produces validated hits in `data/mined_hits/`.

### Option B: Start from pre-mined candidates (skip Phase 0)

Place files directly in `data/mined_hits/`:

<details>
<summary><b>FASTA format</b> — <code>deep_hits_*.fasta</code></summary>

```
>NZ_JAASWF010000012.1_ORF_f1_6252_HIGH
MKISKVDHTRMAVAKGNQHRRDEIGKGLKEVLG...
>NZ_JAAROR010000001.1_ORF_f1_65475_HIGH
MFDKISKVREKNATLKQE...
```

Each sequence is a Cas13 protein ORF. The suffix `_HIGH`/`_MEDIUM`/`_LOW` indicates confidence tier from Phase 0 mining.

</details>

<details>
<summary><b>Metadata CSV format</b> — <code>deep_hits_*_metadata.csv</code></summary>

```csv
sequence_id,repeat_domains,sra_accession,score,confidence,subtype,dr_distance_bp
NZ_JAASWF010000012.1_ORF_f1_6252_HIGH,GTTGTAGCTCCCTTTCTCATTTCGCAGTGCTC,NZ_JAASWF,95.2,HIGH,cas13a,45
NZ_JAPPSR010000004.1_ORF_f1_66985_HIGH,GTCGGCACCGCTCCCGTATAGCGGGG,NZ_JAPPSR,87.5,HIGH,cas13b,112
```

- `repeat_domains`: pipe-separated CRISPR direct repeat k-mers (from adjacent arrays)
- `confidence`: HIGH (BLAST + CRISPR + dual HEPN), MEDIUM (BLAST + dual HEPN), LOW (dual HEPN only)
- `dr_distance_bp`: distance from nearest CRISPR array to the Cas13 gene (38-112 bp is ideal)
- The pipeline selects the optimal repeat per baseline via `fix_crrna_assignments.py`

</details>

---

## 📝 License

This is a research project. Contact the authors for licensing inquiries.

---

<div align="center">

*Built for the frontier of programmable biology. Powered by Rust.*

**CASCADE · *The Dormant Blade Protocol*** — turning nature's collateral curse into medicine's most precise switch.

</div>
