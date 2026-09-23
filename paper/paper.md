---
title: 'CASCADE: oracle-in-the-loop active learning for Cas13 mining and allosteric switch design'
tags:
  - Python
  - CRISPR
  - Cas13
  - protein design
  - structure prediction
  - metagenomics
authors:
  - name: Tyler Waite
    orcid: 0009-0001-0673-2099
    affiliation: "1"
affiliations:
  - name: Independent researcher, San Francisco, CA, USA
    index: 1
date: 22 September 2026
bibliography: paper.bib
---

# Summary

CASCADE (Cas Collateral Activation — Discovery and Engineering) is open-source
software for (1) mining Cas13-like effectors from metagenomic contigs with
strict quality filters and (2) redesigning only inter-domain **linkers** so a
structure model predicts clearer OFF versus ON catalytic geometry. Recognition
(REC / native crRNA) and wild-type HEPN sequences stay frozen. PXDesign and
ProteinMPNN propose linker changes; Protenix or Cattle-Prod score unbound,
trigger-bound, and mismatch complexes; `EvolutionGym` converts relative fitness
into a ProteinMPNN PSSM for the next generation. We demonstrate the stack on
three verified *Listeria booriae* Cas13a proteins that underpin a filed
provisional patent on switch scaffolds, and shortlist linker morphs with
positive composite fitness for continued evolution and wet-lab testing. All
scores here are *in silico*.

# Statement of need

Calling a structure model or inverse folder once is easy; running them as a
**closed, constrained search** is not [@abramson2024af3; @protenix2024;
@dauparas2022mpnn]. Unconstrained redesign often damages catalysis, and
single-state confidence does not measure a switch. CRISPR CAD tools such as
ADAPT [@metsky2022adapt] and BADGERS [@badgers2024] optimize detection guides.
CASCADE targets **effector mechanics**: freeze maps, OFF/ON/mismatch scoring,
hierarchical oracle cost, and a miner that logs rejection reasons so false
Cas13 calls do not poison search [@cascade_audit2026]. Users are computational
biologists and engineers wrapping PXDesign [@pxdesign2025], ProteinMPNN, and
Protenix-compatible oracles (optional Cattle-Prod).

# State of the field

Diffusion binders and global MPNN redesign optimize interfaces or monomer
metrics. CASCADE freezes recognition and catalysis, mutates mechanical
couplers, scores conformational hypotheses, stitches wild-type HEPNs after
generation, and logs JSONL for optional post-training [@wang2024drakes].
Generators and oracles stay external; weights are not redistributed.

# Software design

![CASCADE loop: mine and bootstrap once, then select → generate linkers → stitch HEPNs → hierarchical OFF/ON/mismatch scoring → EvolutionGym PSSM → next generation.](figures/fig_pipeline.png){#fig:loop}

Entrypoint: `scripts/evolution_orchestrator.py` (`scripts/smoke_gpu.sh` for a
short GPU run).

- **Phase 0 — mining (`mining_v3`).** Length gates, canonical `R-X(4-6)-H`
  HEPN pairs, hard anti-signatures, reciprocal Cas13 homology, array-aware
  direct repeats, rejection reason per ORF.
- **Phase 1 — bootstrap.** SQLite annotation, HEPN anchors, native crRNA,
  optional structure screens.
- **Phase 2 — evolution.** Gen 0 scores wild type. Gen 1: PXDesign linker
  backbones. Gen 2+: ProteinMPNN on a cached backbone (≈10–50% linker sites
  free, bias-weighted). Stitch WT HEPNs; mini OFF/ON gates expensive ternary
  and mismatch jobs; relative fitness updates the PSSM.

# Methods

**Mining.** ORFs (default 600–1400 aa) need a valid HEPN pair, no hard
non-Cas13 signatures, and reciprocal hits to `data/cas13_reference_db.fasta`
(remine defaults ≥30% identity / ≥80 aa). CRISPR DRs come from nearby arrays,
not naïve k-mers [@cascade_audit2026].

**Freeze map.** For the three *L. booriae* baselines, HEPN1 is residues
484–594 and HEPN2 is 795–905 (1-based). IDL1 sits before HEPN1; IDL2 between
HEPNs. `hepn_structural_stitch.py` restores WT HEPNs after sampling; stitch
failure costs a fixed fitness penalty (default 5.0).

**Triggers and oracle.** Fusion RNAs in `data/fusion_targets.json` build OFF
(protein+crRNA) and ON (plus target) jobs; 1/2/3-nt spacer mismatches probe
leak. `EVAL_CMD` selects Cattle-Prod or Protenix. His NE2–NE2 distance gates
default to OFF ≥ 18 Å and ON ≤ 12 Å (**design thresholds**, not calibrated
biophysics). Passers may receive base ternary evaluation.

**Generators and gym.** Gen 1: PXDesign. Gen ≥2: ProteinMPNN with unfreeze
`min(0.10+0.05·gen, 0.50)`, `--pssm_multi 0.5`, annealing temperature.
Fitness rewards OFF−ON shift, ipTM, and AF2-IG, doubles for full ternary, and
penalizes mismatch closure. EvolutionGym centers on Gen-0 fitness, EMA-updates
substitution weights (`0.5/0.5`), clips to $[-5,5]$, and exports
`mpnn_bias_gen_*.json`. Elite requires ipTM ≥ 0.85, AF2-IG ≥ 0.80, and ON ≤ 12 Å.
Seeds: `--seed` / `CASCADE_SEED`. CPU tests cover core helpers; figures rebuild
via `scripts/build_paper_figures.py`.

# Demonstration: three Cas13a baselines and morphs

Strict remine kept three *L. booriae* Cas13a ORFs with adjacent CRISPR arrays
[@cascade_audit2026]:

| Baseline (abbrev.) | Contig | aa | Native DR (truncated) | Best ref |
|:-------------------|:-------|---:|:---------------------|:---------|
| NZ_JAASWF…_HIGH | NZ_JAASWF010000012.1 | ~1064 | `UACCUCAAAACAGAAGAGGACUA…` | LseCas13a ~39% |
| NZ_JAAROR…_HIGH | NZ_JAAROR010000001.1 | 1056 | `GAGUACCUCAAAACAGAAGAGGACUAAA…` | LseCas13a ~39% |
| NZ_JAARYF…_HIGH | NZ_JAARYF010000017.1 | 1052 | `GAUUUAGAGUACCUCAAAACAGAAGAGG…` | LseCas13a ~39% |

Campaign 2 (*Bacteroides* / *Flavobacterium* / *Leptotrichia*) accepted
**0/102** ORFs under the same filters (\autoref{fig:mining}). These three
scaffolds are the sequences covered by a **provisional patent application**
filed on the CASCADE switch compositions (sequence listing prepared April
2026); CASCADE is the discovery and design software behind that filing.

![mining_v3 Campaign 1 vs Campaign 2: evaluated vs accepted ORFs, and Campaign 2 rejection reasons.](figures/fig_mining_campaigns.png){#fig:mining}

Intended use (unvalidated): fusion-RNA–gated collateral RNase (e.g. BCR-ABL1,
EWSR1-FLI1) or Cas13 diagnostics. Search concentrated at the HEPN1 border
(**475–484**) and **IDL2 ~700–780** (\autoref{fig:hot}).

![JAAROR: fitness by generation and OFF vs ON distances for top morphs (gates OFF ≥ 18 Å, ON ≤ 12 Å).](figures/fig_evolution_JAAROR_lineage.png){#fig:evo}

![Top-quartile mutation hotspots on JAAROR; blue = IDL2.](figures/fig_mutation_hotspots_JAAROR.png){#fig:hot}

No archived run met elite cuts under the limited GPU budget. Best JAAROR
morph `L59c434_g02_v28` reached fitness ≈ 31.5 (OFF ≈ 66 Å, ON ≈ 20 Å,
ipTM ≈ 0.34); neighbors with large OFF→ON shifts are in
`paper/figures/wetlab_shortlist_JAAROR.csv`. Those candidates—and the three
parental scaffolds—are still being evolved in follow-on CASCADE runs. Absence
of an elite hit reflects incomplete search under funding constraints, not
scaffold exhaustion.

# Research impact statement

CASCADE is already the operational pipeline for the author’s Cas13 switch
program: it produced the three remine-verified *L. booriae* Cas13a baselines
that are the subject of a filed provisional patent on switch scaffolds and
compositions, and it generated the JAAROR linker morphs with the strongest
composite fitness in the archived campaign (e.g. `L59c434_g02_v28` and the
wet-lab shortlist). Those morphs remain under active CASCADE evolution while
GPU time allows. Reproducible mining tables, RL JSONL, figures, Apache-2.0
code, dual-use notes, CPU tests, and the Zenodo software archive
[@cascade_audit2026] document that research use for others to extend.

# AI usage disclosure

Coding assistants (including Cursor) helped with packaging and earlier
manuscript drafts. Methods, figures, and claims were checked by the author
against remine docs and RL logs. Assistants were not used in place of
experiments. Core design decisions (freeze/evolve/score policy, gates, mining
filters) were made by the author.

# Acknowledgements

Built on Protenix, PXDesign, ProteinMPNN, and related open tools. No external
grant funded the evolution campaigns reported here.

# References
