<div align="center">

# 🧬 CASCADE

### **C**as **C**ollateral **A**ctivation — **D**iscovery & **E**ngineering

**Weaponizing Cas13 Collateral Cleavage: An AI-Driven Structural Pipeline for Engineered Suicide Switches in Targeted Oncology**

[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue?logo=python&logoColor=white)](https://www.python.org)
[![Protenix 1.0.4](https://img.shields.io/badge/Protenix-1.0.4-teal?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIyNCIgaGVpZ2h0PSIyNCI+PHBhdGggZD0iTTEyIDJDNi40OCAyIDIgNi40OCAyIDEyczQuNDggMTAgMTAgMTAgMTAtNC40OCAxMC0xMFMxNy41MiAyIDEyIDJ6bTAgMThjLTQuNDIgMC04LTMuNTgtOC04czMuNTgtOCA4LTggOCAzLjU4IDggOC0zLjU4IDgtOCA4eiIgZmlsbD0id2hpdGUiLz48L3N2Zz4=)](https://github.com/bytedance/protenix)
[![PXDesign](https://img.shields.io/badge/PXDesign-0.5.0+pxd-green?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIyNCIgaGVpZ2h0PSIyNCI+PHBhdGggZD0iTTEyIDJDNi40OCAyIDIgNi40OCAyIDEyczQuNDggMTAgMTAgMTAgMTAtNC40OCAxMC0xMFMxNy41MiAyIDEyIDJ6IiBmaWxsPSJ3aGl0ZSIvPjwvc3ZnPg==)](https://github.com/bytedance/PXDesign)
[![Tests](https://img.shields.io/badge/tests-56%20passed-brightgreen?logo=pytest&logoColor=white)](#testing)
[![License](https://img.shields.io/badge/license-Research-lightgrey)](#license)

<br/>

[**Quick Start**](#-quick-start) · [**Architecture**](#-pipeline-architecture) · [**Deploy on VPS**](VPS_DEPLOY.md) · [**RunPod Docker**](RUNPOD_SETUP.md) · [**Interactive Diagram**](workflow_diagram.html) · [**RL Training Data**](RL_TRAINING_FORMAT.md)

<br/>

</div>

---

## 💡 The Core Insight

> Traditional gene therapies view Cas13's **collateral trans-cleavage** — the indiscriminate shredding of bystander RNA upon activation — as a catastrophic flaw. **CASCADE reverses this paradigm.**

We discover, validate, and engineer novel Cas13e-like proteins from metagenomic dark matter to function as **highly specific biological suicide switches**:

| State | HEPN Domains | Behavior |
|:------|:-------------|:---------|
| 🟢 **Healthy cell** (OFF) | **> 25 Å apart** — catalytically inert | Completely dormant. Zero leakiness. |
| 🔴 **Tumor cell** (ON) | **< 12 Å apart** — catalytically aligned | Massive collateral cleavage → programmed cell death |

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

CASCADE operates in two phases. Phase 1 bootstraps validated baselines from raw metagenomic data. Phase 2 runs an autonomous **reinforcement-learning (RL) driven evolution loop** that continuously designs, evaluates, and learns from structural predictions.

```
  ┌─────────────────────────────────────────────────────────────┐
  │  PHASE 1: Bootstrap                                         │
  │  Raw FASTAs → Parse & HEPN-anchor → Protenix-mini screen   │
  │  → Validated baselines with native crRNAs                   │
  └────────────────────────┬────────────────────────────────────┘
                           ▼
  ┌─────────────────────────────────────────────────────────────┐
  │  PHASE 2: Active Learning Evolution Loop                    │
  │                                                             │
  │  ┌──────────┐    ┌──────────┐    ┌──────────┐              │
  │  │ PXDesign │───▶│ Protenix │───▶│ Fitness  │              │
  │  │ Generate │    │ Evaluate │    │  Score   │              │
  │  └────▲─────┘    └──────────┘    └────┬─────┘              │
  │       │                               │                     │
  │       │    ┌──────────────────┐       │                     │
  │       └────│  RL Bias Matrix  │◀──────┘                     │
  │            │  (EvolutionGym)  │                              │
  │            └──────────────────┘                              │
  │                                                             │
  │  Global Best Tracking: evolution ALWAYS proceeds from the   │
  │  highest-fitness protein discovered so far — never regresses│
  └─────────────────────────────────────────────────────────────┘
```

> 📊 **[Open the interactive workflow diagram →](workflow_diagram.html)** for the full node-by-node breakdown with color-coded phases.

### Phase 1 — Bootstrap Baselines

| Step | Script | What It Does |
|:-----|:-------|:-------------|
| **1a** | `01_parse_and_annotate.py` | Parse FASTA + CSV into SQLite DB. Anchor HEPN1/HEPN2 domains via `R.{4,6}H` motif. Generate Protenix-compatible JSONs. |
| **1b** | `02_run_screening.sh` | GPU-accelerated Protenix-mini structural screen. Filters hits that can't form bilobed structures or bind crRNA. |
| **1c** | `02b_rerun_top_with_msa.sh` | *(Optional)* Re-run top N baselines with MSA for higher-quality seed structures. |
| **1d** | `validate_crispr_repeats.py` | *(Optional)* Validate CRISPR repeats via RNAfold. Outputs `validated_baseline_ids.txt` to restrict evolution. |

### Phase 2 — Active Learning Evolution Loop

The evolution orchestrator (`evolution_orchestrator.py`) runs an autonomous loop:

```
For each generation:
  1. GENERATE  → PXDesign designs linker variants (REC + HEPN frozen)
  2. STITCH    → Wild-type HEPN1/HEPN2 grafted into designed linkers
  3. RESOLVE   → Unknown residues (X) replaced from baseline or Glycine
  4. BIAS      → RL bias matrix applied to linker regions (closed-loop RL)
  5. EVALUATE  → Protenix predicts OFF-state and ON-state structures
  6. MEASURE   → 3D HEPN distance (Å) + ipTM + AF2-IG confidence scores
  7. TEST      → Off-target specificity (1/2/3-mismatch progressive penalty)
  8. SCORE     → Composite fitness: HEPN_shift + ipTM×50 + AF2-IG×20 − penalties
  9. LEARN     → EvolutionGym updates mutation weights → new bias matrix
  10. ADVANCE  → Global best protein becomes next baseline (never regresses)
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

**Thresholds:**
- OFF distance ≥ 25 Å (catalytically inert when dormant)
- ON distance ≤ 12 Å (catalytically aligned when activated)
- Elite: ipTM ≥ 0.85, AF2-IG ≥ 0.80, ON ≤ 12 Å

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

## 🏷️ Variant Naming Convention

Every variant produced by CASCADE follows a fixed naming scheme designed to stay unique and sortable across thousands of runs:

```
R{RRR}.L{LL}.G{GG}.V{VVV}
```

| Segment | Meaning | Padding | Example |
|:--------|:--------|:--------|:--------|
| `R` | **Run** — independent pipeline execution | 3 digits | `R001` |
| `L` | **Lineage** — baseline lineage within a run | 2 digits | `L03` |
| `G` | **Generation** — evolution generation within a lineage | 2 digits | `G12` |
| `V` | **Variant** — design index within a generation | 3 digits | `V007` |

**Suffix modifiers** (appended after the variant segment):

| Suffix | Meaning |
|:-------|:--------|
| `F` | Fallback variant (stitch failure → baseline reused, penalty applied) |
| `E` | Elite variant (passed all elite thresholds) |

**Examples:**

| Name | Meaning |
|:-----|:--------|
| `R001.L02.G05.V023` | Run 1, lineage 2, generation 5, variant 23 |
| `R001.L02.G05.V023F` | Same as above — fallback (stitching failed) |
| `R014.L01.G19.V004E` | Run 14, lineage 1, generation 19, variant 4 — elite |

> Zero-padded segments keep filenames lexicographically sorted. The scheme supports up to 999 runs × 99 lineages × 99 generations × 999 variants per generation.

---

## 🚀 Quick Start

### Prerequisites

| Requirement | Minimum |
|:------------|:--------|
| **GPU** | A100 40GB (80GB recommended) |
| **OS** | Ubuntu 22.04+ with CUDA 12.x |
| **Python** | 3.11+ |
| **Conda** | Miniconda / Anaconda |
| **Storage** | 100–200 GB free |

### Setup (One-Time)

```bash
git clone https://github.com/twaite11/CASCADE-Cas-Collateral-Activation-Discovery-Engineering.git
cd CASCADE-Cas-Collateral-Activation-Discovery-Engineering

# Automated dual-environment setup (Protenix 1.0.4 + PXDesign/Protenix 0.5.0)
chmod +x scripts/setup_dual_env.sh
./scripts/setup_dual_env.sh
```

This creates two isolated conda environments:

| Environment | Protenix | Purpose |
|:------------|:---------|:--------|
| `cascade` | 1.0.4 | Structure evaluation, fitness scoring, evolution loop |
| `pxdesign` | 0.5.0+pxd | PXDesign variant generation (called cross-env) |

> **Why two envs?** PXDesign depends on Protenix 0.5.0+pxd, which conflicts with Protenix 1.0.4. The pipeline runs in `cascade` and calls PXDesign from `pxdesign` via `PXDESIGN_CMD`.

### Activate (Every Session)

```bash
source scripts/cascade_env.sh
```

### Run

```bash
# Place your data in data/mined_hits/
#   deep_hits_*.fasta     — Cas13e-like protein sequences
#   deep_hits_*_metadata.csv — Metadata with crRNA repeat domains

cd scripts

# Phase 1: Parse → Screen
python 01_parse_and_annotate.py
./02_run_screening.sh

# Phase 2: Evolution
python evolution_orchestrator.py
```

Or run everything at once:

```bash
./scripts/run_pipeline.sh 2>&1 | tee "logs/cascade_$(date +%Y%m%d_%H%M%S).log"
```

> 📖 **Full deployment guide:** [VPS_DEPLOY.md](VPS_DEPLOY.md)

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
│
├── 📁 data/
│   └── 📁 mined_hits/                  # Input: FASTAs + metadata CSVs
│       ├── deep_hits_*.fasta           #   Cas13e-like protein ORFs
│       └── deep_hits_*_metadata.csv    #   SRA accessions + crRNA k-mers
│
├── 📁 scripts/
│   ├── 🔧 setup_dual_env.sh           # Creates cascade + pxdesign conda envs
│   ├── 🔧 setup_vps.sh                # Simple single-env venv setup
│   ├── 🔧 cascade_env.sh              # Source this to activate env + PXDESIGN_CMD
│   ├── 🔧 run_pipeline.sh             # Full pipeline runner with logging
│   │
│   ├── 📜 01_parse_and_annotate.py     # Phase 1a: Ingest → SQLite DB + HEPN anchoring
│   ├── 📜 02_run_screening.sh          # Phase 1b: Protenix-mini structural screen
│   ├── 📜 02b_rerun_top_with_msa.sh   # Phase 1c: Optional MSA re-run for top N
│   ├── 📜 validate_crispr_repeats.py   # Phase 1d: Optional CRISPR repeat validation
│   │
│   ├── 📜 03_pxdesign_wrapper.py       # Variant generation (PXDesign + stitching + RL bias)
│   ├── 📜 evolution_orchestrator.py    # Phase 2: Active learning master controller
│   │
│   └── 📁 utils/
│       ├── 📜 protenix_eval.py         # ON/OFF payload generation + Protenix inference
│       ├── 📜 pdb_kinematics.py        # 3D HEPN distance + confidence score extraction
│       └── 📜 hepn_structural_stitch.py # Graft WT HEPN domains into designed linkers
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
│   ├── cas13_variants.db               #   SQLite database
│   └── variant_domain_metadata.json    #   HEPN domain boundaries
│
├── 📁 jsons/                           # Protenix-compatible input payloads
│
└── 📁 outputs/
    ├── phase1_screening/               # Baseline CIF/PDB structures
    ├── generation_queue/               # PXDesign variant FASTAs
    ├── fast_eval/                      # Mini-model OFF/ON screening
    ├── high_fidelity_scoring/          # Base-model ternary predictions
    ├── rl_gym_data/                    # RL bias matrices + training dataset
    │   ├── mpnn_bias_gen_X.json        #   Per-generation bias matrices
    │   └── rl_training_dataset.jsonl   #   For DRAKES/ProteinMPNN post-training
    └── optimized_switches/             # 🏆 Elite outputs
        ├── *_optimal.fasta             #   Best protein sequences
        ├── *_ternary_complex.cif       #   Predicted ternary structures
        └── *_crRNA.fasta               #   Native crRNA sequences
```

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
| `MAX_GENERATIONS` | 20 | Generations per lineage |
| `NUM_INITIAL_LINEAGES` | 5 | Parallel lineages from top Phase 1 baselines |
| `MIN_OFF_DISTANCE` | 25.0 Å | Minimum HEPN distance for dormant state |
| `MAX_ON_DISTANCE` | 12.0 Å | Maximum HEPN distance for active state |
| `FALLBACK_FITNESS_PENALTY` | 5.0 | Penalty for stitch-failure fallback variants |
| `SLEEP_AFTER_PXDESIGN` | 5 s | GPU cooldown between PXDesign and Protenix |
| `SLEEP_AFTER_PROTENIX` | 2 s | GPU cooldown between Protenix runs |

Environment variables:

| Variable | Example | Description |
|:---------|:--------|:------------|
| `PXDESIGN_CMD` | `conda run --no-banner -n pxdesign pxdesign` | How to invoke PXDesign from the cascade env |
| `PROTENIX_BASE_MODEL` | `protenix_base_default_v1.0.0` | Override Protenix base model name |
| `CUDA_VERSION` | `12.1` | CUDA version for dual-env setup |

---

## 🧬 Expected Input Data

Place files in `data/mined_hits/`:

<details>
<summary><b>FASTA format</b> — <code>deep_hits_*.fasta</code></summary>

```
>Cas13a_positive_control
MKISKVDHTRMAVAKGNQHRRDEIGKGLKEVLG...
>SRA_hit_001
MFDKISKVREKNATLKQE...
```

Each sequence is a predicted Cas13e-like protein ORF from metagenomic contigs.

</details>

<details>
<summary><b>Metadata CSV format</b> — <code>deep_hits_*_metadata.csv</code></summary>

```csv
sequence_id,repeat_domains,sra_accession,score
Cas13a_positive_control,GATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC,known,95.2
SRA_hit_001,GTTGTAGCTCCCTTTCTCATTTCGCAGTGCTC|GTTGTAGCTCCCTTACTCATTTCGGAGTGCTC,SRR12345678,87.5
```

- `repeat_domains`: pipe-separated CRISPR direct repeat k-mers
- The pipeline selects the optimal repeat per baseline

</details>

---

## 📝 License

This is a research project. Contact the authors for licensing inquiries.

---

<div align="center">

*Built for the frontier of programmable biology.*

**CASCADE** — turning nature's "flaw" into medicine's most precise weapon.

</div>
