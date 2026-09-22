# CASCADE

**Cas** **C**ollateral **A**ctivation — **D**iscovery & **E**ngineering

An open computational pipeline for (1) strict Cas13-like effector mining from
metagenomic contigs, (2) baseline bootstrap with native CRISPR repeats and
structure screens, and (3) **oracle-in-the-loop** constrained sequence search
over inter-domain linkers, scored by predicted OFF vs ON catalytic geometry.

[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22905895.svg)](https://doi.org/10.5281/zenodo.22905895)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue?logo=python&logoColor=white)](https://www.python.org)
[![JOSS](https://img.shields.io/badge/JOSS-submitted-yellow)](paper/paper.md)

[Quick start](#quick-start) · [GPU smoke](#gpu-smoke-recipe) · [Architecture](#pipeline-architecture) · [Dual-use](docs/DUAL_USE.md) · [Third-party](THIRD_PARTY_NOTICES.md) · [JOSS paper](paper/paper.md) · [RL data format](RL_TRAINING_FORMAT.md)

---

## What CASCADE is (methods)

CASCADE freezes RNA recognition (REC / native crRNA) and catalytic HEPN
domains, then searches **mechanical** regions (inter-domain linkers) with
PXDesign / ProteinMPNN. A hierarchical structure oracle (Cattle-Prod or
Protenix) evaluates:

| State | Job | Design proxy |
|:------|:----|:-------------|
| **OFF** | Protein + crRNA, no cognate trigger | HEPN His NE2–NE2 distance high (default ≥ 18 Å) |
| **ON** | Protein + crRNA + trigger RNA | Distance low (default ≤ 12 Å), plus ipTM / AF2-IG |
| **Mismatch** | ON jobs with 1/2/3-nt spacer mismatches | Penalize activation on near-cognate RNA |

An online adapter (`EvolutionGym`) converts relative fitness into a ProteinMPNN
PSSM bias for the next generation and writes a `(sequence, structure, reward)`
JSONL dataset for optional offline fine-tuning.

**Evidence class today:** *in silico* only. Distance gates are design
thresholds, not wet-lab calibrated catalytic measurements. See
[Intended application (unvalidated)](#intended-application-unvalidated) and
[docs/DUAL_USE.md](docs/DUAL_USE.md).

---

## Intended application (unvalidated)

> This section is a **research hypothesis**, not a validated result.

Cas13 collateral *trans*-cleavage is usually treated as a diagnostic side
effect. One possible application of a gated effector is a
**tumor-restricted RNA-triggered switch**: activation only when a fusion RNA
(e.g. BCR-ABL1, EWSR1-FLI1 from `data/fusion_targets.json`) is bound, with
HEPN domains predicted to assemble into a cleavage-competent geometry inside
that cell type and remain apart otherwise.

CASCADE does **not** claim:

- novel “Cas13e-like dark matter” therapeutics (the audited catalog is three
  *Listeria booriae* Cas13a homologs — see `docs/REMINE_2026-05-13.md`)
- zero leakiness, programmed cell death, or clinical utility
- freedom-to-operate against patented Cas13 sequences

Wet-lab validation (expression, collateral RNase assays, ON vs OFF transcripts)
is out of scope for this software release.

---

## Engineering contract (freeze / evolve / score)

| Do not touch | Search this |
|:-------------|:------------|
| REC / crRNA pocket — use the **native** direct repeat from the contig | Inter-domain linkers (and related mechanical segments) |
| Wild-type HEPN1 / HEPN2 catalytic sequences (stitched back after generation) | Sequence / backbone proposals from PXDesign (Gen 1) and ProteinMPNN (Gen 2+) |

Fork surface for other enzymes: replace freeze map, catalytic atom pairs,
trigger job templates, and fitness weights (today these are still
Cas13-shaped; see roadmap in `paper/paper.md`).

---

## Quick start

### Prerequisites

| Requirement | Minimum |
|:------------|:--------|
| GPU (Phase 1–2) | A100 40GB recommended |
| OS | Ubuntu 22.04+ with CUDA 12.x (CPU tests need no GPU) |
| Python | 3.11+ |
| Rust (optional) | 1.78+ for Cattle-Prod / `rust/` accelerators |

**Weights are not in this repo.** Install Protenix and/or Cattle-Prod, PXDesign,
and ProteinMPNN yourself. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

### Setup

```bash
git clone https://github.com/twaite11/CASCADE-Cas-Collateral-Activation-Discovery-Engineering.git
cd CASCADE-Cas-Collateral-Activation-Discovery-Engineering

# Optional: Cattle-Prod eval engine
git clone https://github.com/twaite11/cattle-prod.git ../cattle-prod
# build per cattle-prod docs, then:
export PATH="$(realpath ../cattle-prod/cattle-prod/target/release):$PATH"

chmod +x scripts/setup_dual_env.sh
./scripts/setup_dual_env.sh
source scripts/cascade_env.sh
```

### CPU tests (no GPU)

```bash
pip install -r requirements.txt pytest
pytest tests/ -q
```

### Phase 0 — mining (CPU-friendly)

```bash
# Place contigs in data/new_contigs/, then:
python scripts/mining_v3.py \
  --contigs data/new_contigs/*.fasta \
  --output-dir outputs/mining_v3
```

Confirmed baselines used in docs live under `data/mined_hits/confirmed_novel_cas13.*`.

### Phase 1–2 — science path (local GPU, no Vast.ai)

```bash
# Bootstrap annotations + structural screen (see scripts/02_*.sh)
python scripts/01_parse_and_annotate.py
# ... then screening scripts once EVAL_CMD is set ...

# Evolution loop (primary scientific entrypoint)
CASCADE_WORKERS=1 python scripts/evolution_orchestrator.py \
  --baseline-ids NZ_JAASWF010000012.1_ORF_f1_6252_HIGH \
  --max-generations 12 \
  --seed 42
```

Or use the short smoke recipe:

```bash
chmod +x scripts/smoke_gpu.sh
./scripts/smoke_gpu.sh
```

### Optional: dashboard + Vast.ai

Cloud fan-out is **optional ops**, not required for science. See
[CONTAINER_DEPLOY.md](CONTAINER_DEPLOY.md) and [VPS_DEPLOY.md](VPS_DEPLOY.md).
Local science should always work via `evolution_orchestrator.py` +
`scripts/smoke_gpu.sh`.

---

## Pipeline architecture

```
Phase 0  mining_v3          → ranked Cas13-like ORFs + DRs + rejection_reason
Phase 1  annotate + screen  → SQLite, HEPN anchors, mini/base structures
Phase 2  evolution loop     → generate → stitch → OFF/ON/mismatch eval → gym PSSM
```

Interactive diagram: [workflow_diagram.html](workflow_diagram.html).

### Evolution loop (AI engineering)

1. **Gen 0** — score unmodified enzyme (fitness baseline)
2. **Gen 1** — PXDesign proposes linker geometry
3. **Stitch** — graft WT HEPN domains; fail closed with penalty
4. **Mini OFF/ON** — cheap oracle; gate on distance thresholds
5. **Base ternary + mismatches** — expensive oracle for passers
6. **EvolutionGym** — relative advantage → EMA mutation weights → MPNN PSSM
7. **Gen 2+** — ProteinMPNN on cached backbone; unfreeze 10–50% of linker sites
8. **JSONL flywheel** — append reward-labeled rows for offline fine-tuning

Fitness (simplified): HEPN conformational shift + confidence terms − mismatch
penalties. Constants live in `scripts/evolution_orchestrator.py` and are
overridable via environment variables where noted.

---

## License and third-party stack

- **CASCADE code:** [Apache License 2.0](LICENSE) — see also [NOTICE](NOTICE)
- **Runtime dependencies / weights:** [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)

  | Tool | Code license | Weights |
  |:-----|:-------------|:--------|
  | Protenix | Apache-2.0 | Not redistributed — download upstream |
  | PXDesign | Apache-2.0 | Not redistributed — download upstream |
  | Cattle-Prod | Apache-2.0 | Not redistributed — convert/place yourself |
  | ProteinMPNN | MIT | Not redistributed — use upstream |

---

## Frozen results (Zenodo)

Mining tables and curated GPU artifacts should be archived on Zenodo (not only
under gitignored `outputs/`).

```bash
python scripts/build_zenodo_bundle.py
# optional larger zip with one CIF sample per baseline:
python scripts/build_zenodo_bundle.py --include-cif-samples
```

Upload `dist/zenodo/CASCADE_results_*.zip` at
[zenodo.org/deposit/new](https://zenodo.org/deposit/new) (CC-BY-4.0 for data).
Then paste the DOI into `CITATION.cff` and `paper/paper.md`.

Committed reproducible summaries also live in `outputs/mining_v3_campaign*/`
and `docs/REMINE_2026-05-13.md`.

---

## Dual-use

Tumor-restricted fusion RNA targeting is the **intended research framing**.
All scoring is **in silico**. CASCADE is **not a therapeutic**. Full statement:
[docs/DUAL_USE.md](docs/DUAL_USE.md).

---

## Citation

See [CITATION.cff](CITATION.cff). JOSS draft: [paper/paper.md](paper/paper.md).

---

## Project layout (short)

```
scripts/mining_v3.py              # Strict Cas13 mining
scripts/evolution_orchestrator.py # Primary evolution entrypoint
scripts/smoke_gpu.sh              # One-baseline GPU smoke
scripts/build_zenodo_bundle.py    # Results archive for Zenodo
scripts/utils/                    # Oracle I/O, HEPN stitch, kinematics
paper/                            # JOSS manuscript
dashboard_backend/                # Optional FastAPI UI
tests/                            # CPU unit tests
```

---

*CASCADE — constrained search around frozen structure oracles. Computational methods software; application claims unvalidated.*
