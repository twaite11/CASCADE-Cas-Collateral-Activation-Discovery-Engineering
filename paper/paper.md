---
title: 'CASCADE: oracle-in-the-loop active learning for constrained Cas13 switch design'
tags:
  - Python
  - Rust
  - protein design
  - CRISPR
  - Cas13
  - structure prediction
  - active learning
authors:
  - name: Tyler Waite
    # orcid: 0000-0000-0000-0000   # uncomment with your ORCID before JOSS submit
    affiliation: "1"
affiliations:
  - name: Independent researcher, San Francisco, CA, USA
    index: 1
date: 22 September 2026
bibliography: paper.bib
---

# Summary

CASCADE (Cas Collateral Activation — Discovery and Engineering) is open-source
software for discovering Cas13-like effectors from metagenomic contigs and for
running a constrained, oracle-in-the-loop sequence search over inter-domain
linkers. Recognition and catalytic domains are frozen; generative models
propose only mechanical segments. A hierarchical structure-prediction oracle
scores unbound (OFF) versus trigger-bound (ON) complexes and near-cognate
mismatch controls. An online adapter converts relative fitness into a
ProteinMPNN position-specific scoring matrix for the next generation and
writes a reward-labeled dataset for optional offline fine-tuning. The
reference instance uses tumor fusion RNA triggers as an *in silico*
hypothesis; CASCADE does not claim wet-lab or clinical validation.

# Statement of need

Foundation models for biomolecular structure (@abramson2024af3; @protenix2024)
and inverse folding (@dauparas2022mpnn) are widely available, but turning them
into a reliable design loop requires systems engineering: constrained action
spaces so catalytic sites are not destroyed, multi-state evaluation so
“foldedness” is not mistaken for switching, cost-aware cascades so expensive
ternary predictions are reserved for promising candidates, and feedback so
search improves within a run. Existing CRISPR CAD tools such as ADAPT
(@metsky2022adapt) and BADGERS (@badgers2024) optimize guide sequences for
detection. CASCADE addresses a complementary problem: **constrained redesign
of effector mechanics** under two-state structural scoring, plus a strict
mining stage that documents rejection reasons when putative Cas13 hits are
false positives (@cascade_audit2026).

Target users are computational biologists and AI engineers who need a
reproducible pipeline—not only notebooks—around PXDesign (@pxdesign2025),
ProteinMPNN, and Protenix-compatible oracles (including the optional Rust
Cattle-Prod backend).

# State of the field

Unconstrained protein design stacks (diffusion binders, global MPNN redesign)
optimize for interface confidence or monomer metrics. CASCADE instead encodes
an explicit **freeze / evolve / score** policy for multi-domain, ligand-gated
enzymes. Relative to CRISPR diagnostic design packages, CASCADE does not
replace guide activity models; it wraps structure oracles and generators to
search linker sequence space under OFF/ON/mismatch objectives. Relative to
generic active-learning wrappers, CASCADE ships domain-aware stitching of
wild-type HEPN segments, Cas13 mining with anti-signature filters, and a
JSONL flywheel aligned with discrete-diffusion post-training formats
(@wang2024drakes).

# Software design

CASCADE separates concerns so mining quality, oracle cost, and search feedback
can evolve independently:

1. **Mining (`mining_v3`)** — canonical HEPN motifs, hard disqualifiers
   (e.g. GH3, MetH, AAA+), reciprocal homology, CRISPR array–aware direct
   repeat assignment, and per-ORF rejection logs. Earlier looser mining
   produced large false-positive catalogs; the strict path is therefore part
   of the machine-learning story (garbage-in destroys any oracle loop).
2. **Bootstrap** — annotation into SQLite, HEPN anchoring, optional mini/base
   structural screens that produce the starting population for evolution.
3. **Evolution orchestrator** — Gen-0 baseline scoring; Gen-1 PXDesign
   backbones; Gen-2+ ProteinMPNN on a cached backbone with a scheduled
   unfreeze of linker positions (10–50\%); HEPN stitch that fails closed with
   an explicit fitness penalty; hierarchical OFF/ON then ternary/mismatch
   evaluation via `EVAL_CMD` (Cattle-Prod or Protenix).
4. **EvolutionGym** — relative advantage versus the Gen-0 baseline (or
   generation mean), EMA-updated mutation weights, clipped PSSM export into
   ProteinMPNN (`--pssm_multi`). This is online combinatorial adaptation, not
   a learned policy-gradient agent; naming it clearly matters for AI-facing
   readers.
5. **Optional ops** — FastAPI dashboard and Vast.ai fan-out for parallel
   hosts. Scientific runs are intended to work from
   `scripts/evolution_orchestrator.py` and `scripts/smoke_gpu.sh` on a single
   GPU without cloud orchestration.

Trade-offs: CASCADE reuses frozen foundation models as judges rather than
training a new generator; fitness coefficients and Ångström gates are
explicit design thresholds and must not be read as calibrated biophysics;
topology is still Cas13-shaped (REC–linker–HEPN–linker–HEPN), with a
documented path toward a declarative switch specification so other
ligand-gated enzymes can reuse the same freeze/evolve/score contract.

# Research impact statement

The software includes an extensive CPU unit-test suite covering fitness,
stitching, kinematics payloads, determinism hooks, and dashboard/ops helpers;
dual-environment setup scripts; Apache-2.0 licensing with explicit third-party
weight notices; and dual-use documentation stating tumor-restricted RNA
targeting as an *in silico* research framing only. Frozen mining campaign
summaries in-repo (`outputs/mining_v3_campaign*`, `docs/REMINE_2026-05-13.md`)
document recovery of three verified *Listeria booriae* Cas13a effectors and
rejection of 102/102 Campaign-2 ORFs under strict filters—a reproducible
negative control. A builder script (`scripts/build_zenodo_bundle.py`) packages
mining tables and curated GPU artifacts (reward JSONL, bias matrices,
optimized FASTA, Phase-1 confidence summaries) for Zenodo deposit; the DOI
should be linked in `CITATION.cff` after upload. Together these materials
support JOSS-style review: installable code, tests, example baselines, and
archived results separate from foundation-model weights.

# AI usage disclosure

Generative AI coding assistants (including Cursor) were used to help draft
documentation, packaging scripts, and this manuscript text. Pipeline logic,
thresholds, and scientific claims were reviewed and edited by the author.
CPU tests and mining audits described in the repository documentation were
used to check correctness of software behavior. No generative model was used
as an unsupervised substitute for experimental validation.

# Acknowledgements

CASCADE builds on open Protenix, PXDesign, ProteinMPNN, and related
ecosystem tools. No specific funding grant is claimed for this software
release.

# References
