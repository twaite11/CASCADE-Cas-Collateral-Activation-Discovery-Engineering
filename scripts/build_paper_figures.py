#!/usr/bin/env python3
"""Build paper figures from mining + RL campaign data."""
from __future__ import annotations

import json
import shutil
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# Desktop patent figures (if present)
DESK_FIGS = Path.home() / "Desktop" / "patent_figures"


def load_jsonl(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            rows.append(json.loads(line))
    return rows


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.grid(axis="y", linestyle=":", alpha=0.4)


def fig_mining():
    """Campaign 1 vs 2 acceptance (mining_v3)."""
    c1 = json.loads((ROOT / "outputs/mining_v3_campaign1/mining_v3_summary.json").read_text())
    c2 = json.loads((ROOT / "outputs/mining_v3_campaign2/mining_v3_summary.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4), constrained_layout=True)

    ax = axes[0]
    cats = ["Campaign 1\nL. booriae", "Campaign 2\nBacteroides /\nFlavo / Lepto"]
    evaluated = [c1["orfs_evaluated"], c2["orfs_evaluated"]]
    accepted = [c1["accepted"], c2["accepted"]]
    x = np.arange(len(cats))
    ax.bar(x - 0.18, evaluated, width=0.36, label="ORFs evaluated", color="#6b7280")
    ax.bar(x + 0.18, accepted, width=0.36, label="Accepted Cas13", color="#2563eb")
    ax.set_xticks(x)
    ax.set_xticklabels(cats, fontsize=8)
    ax.set_ylabel("Count")
    ax.set_title("Strict mining_v3: evaluated vs accepted")
    ax.legend(fontsize=8, frameon=False)
    style_ax(ax)

    ax = axes[1]
    reasons = c2.get("rejection_counts", {})
    labels = [
        "No R-X(4-6)-H",
        "No reciprocal hit",
        "GH3",
        "Cas9 RuvC",
        "MBL/RNase Z",
    ]
    keys = [
        "no_canonical_R-X(4-6)-H_pair",
        "no_reciprocal_hit(>=30%id over 80aa)",
        "hard_signature:GH3_KHFPGHGD,GH3_CATALYTIC",
        "hard_signature:CAS9_RUVC",
        "hard_signature:MBL_RNASE_Z",
    ]
    vals = [reasons.get(k, 0) for k in keys]
    y = np.arange(len(labels))
    ax.barh(y, vals, color="#b45309")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Rejected ORFs")
    ax.set_title("Campaign 2 rejection reasons (n=102)")
    style_ax(ax)

    out = FIG / "fig_mining_campaigns.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def fig_pipeline():
    """Simple pipeline diagram as matplotlib boxes."""
    fig, ax = plt.subplots(figsize=(8.5, 2.6), constrained_layout=True)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3)
    ax.axis("off")
    boxes = [
        (0.3, 1.0, "Mine\ncontigs", "#e5e7eb"),
        (2.3, 1.0, "Bootstrap\nbaselines", "#e5e7eb"),
        (4.3, 1.0, "Generate\nlinkers only", "#dbeafe"),
        (6.3, 1.0, "Score\nOFF / ON /\nmismatch", "#dbeafe"),
        (8.3, 1.0, "Update\nPSSM bias", "#dbeafe"),
    ]
    for x, y, text, color in boxes:
        ax.add_patch(
            plt.Rectangle((x, y), 1.5, 1.2, facecolor=color, edgecolor="#111", linewidth=1.2)
        )
        ax.text(x + 0.75, y + 0.6, text, ha="center", va="center", fontsize=8)
    for i in range(4):
        x0 = boxes[i][0] + 1.5
        x1 = boxes[i + 1][0]
        ax.annotate(
            "",
            xy=(x1, 1.6),
            xytext=(x0, 1.6),
            arrowprops=dict(arrowstyle="->", color="#111", lw=1.2),
        )
    ax.text(
        5,
        0.35,
        "Freeze REC + HEPNs · mutate linkers · structure oracle is the judge",
        ha="center",
        fontsize=9,
        color="#374151",
    )
    ax.set_title("CASCADE loop (single GPU; Vast.ai optional)", fontsize=11, pad=8)
    out = FIG / "fig_pipeline.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def fig_evolution(rows: list[dict], tag: str):
    """Fitness over generations + OFF/ON scatter for top morphs."""
    by_gen = defaultdict(list)
    for r in rows:
        by_gen[r["generation"]].append(r["fitness"])
    gens = sorted(by_gen)
    med = [statistics.median(by_gen[g]) for g in gens]
    mx = [max(by_gen[g]) for g in gens]

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5), constrained_layout=True)
    ax = axes[0]
    ax.plot(gens, mx, "o-", color="#2563eb", label="Best fitness / gen")
    ax.plot(gens, med, "s--", color="#6b7280", label="Median fitness / gen")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Composite fitness")
    ax.set_title(f"Evolution progress ({tag})")
    ax.legend(fontsize=8, frameon=False)
    style_ax(ax)

    ax = axes[1]
    top = sorted(rows, key=lambda r: r["fitness"], reverse=True)[:40]
    offs = [r["off_dist_A"] for r in top]
    ons = [r["on_dist_A"] for r in top]
    fits = [r["fitness"] for r in top]
    sc = ax.scatter(offs, ons, c=fits, cmap="viridis", s=36, edgecolors="#111", linewidths=0.3)
    ax.axvline(18, color="#9ca3af", ls=":", lw=1)
    ax.axhline(12, color="#9ca3af", ls=":", lw=1)
    ax.set_xlabel("OFF HEPN distance (Å)")
    ax.set_ylabel("ON HEPN distance (Å)")
    ax.set_title("Top morphs: OFF vs ON distance")
    # annotate best
    best = top[0]
    ax.annotate(
        best["variant_id"],
        (best["off_dist_A"], best["on_dist_A"]),
        textcoords="offset points",
        xytext=(6, 6),
        fontsize=7,
    )
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Fitness", fontsize=8)
    style_ax(ax)

    out = FIG / f"fig_evolution_{tag}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out, best


def fig_mutation_map(rows: list[dict], coords: dict | None, tag: str):
    """Hotspot histogram along sequence for top-quartile fitness variants."""
    fits = sorted(r["fitness"] for r in rows)
    thr = fits[int(0.75 * (len(fits) - 1))]
    hot = Counter()
    for r in rows:
        if r["fitness"] < thr:
            continue
        for m in r.get("mutations") or []:
            parts = m.split("_")
            if len(parts) == 2 and parts[1] != "del" and parts[0].isdigit():
                hot[int(parts[0])] += 1
    if not hot:
        return None
    positions = sorted(hot)
    counts = [hot[p] for p in positions]

    fig, ax = plt.subplots(figsize=(8.2, 2.8), constrained_layout=True)
    ax.bar(positions, counts, width=1.0, color="#2563eb", alpha=0.85)
    if coords:
        # Accept either flat coords or nested domains{} from variant_domain_metadata.json
        domains = coords.get("domains", coords)
        h1s = int(domains.get("HEPN1", {}).get("start", coords.get("hepn1_start", 0)) or 0)
        h1e = int(domains.get("HEPN1", {}).get("end", coords.get("hepn1_end", 0)) or 0)
        h2s = int(domains.get("HEPN2", {}).get("start", coords.get("hepn2_start", 0)) or 0)
        h2e = int(domains.get("HEPN2", {}).get("end", coords.get("hepn2_end", 0)) or 0)
        # N-terminal block (REC + IDL1) ends at HEPN1 start; IDL2 sits between HEPNs.
        bands = [
            ("REC+IDL1", 0, h1s, "#e5e7eb"),
            ("HEPN1", h1s, h1e, "#f3f4f6"),
            ("IDL2", h1e, h2s, "#bfdbfe"),
            ("HEPN2", h2s, h2e, "#f3f4f6"),
        ]
        ymax = max(counts) * 1.15
        for name, a, b, color in bands:
            if b > a:
                ax.axvspan(a, b, color=color, alpha=0.55, zorder=0)
                ax.text((a + b) / 2, ymax * 0.92, name, ha="center", fontsize=7, color="#374151")
        ax.set_ylim(0, ymax)
    ax.set_xlabel("Residue position (1-based)")
    ax.set_ylabel("Count in top-quartile variants")
    ax.set_title(f"Mutation hotspots under linker-only search ({tag})")
    style_ax(ax)
    out = FIG / f"fig_mutation_hotspots_{tag}.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def copy_patent_svgs():
    copied = []
    if not DESK_FIGS.is_dir():
        return copied
    mapping = {
        "Figure1_CRISPR_Loci.svg": "fig_crispr_loci.svg",
        "Figure2_HEPN_Architecture.svg": "fig_hepn_architecture.svg",
    }
    for src_name, dest_name in mapping.items():
        src = DESK_FIGS / src_name
        if src.exists():
            dest = FIG / dest_name
            shutil.copy2(src, dest)
            copied.append(dest)
    return copied


def main():
    print("mining", fig_mining())
    print("pipeline", fig_pipeline())
    copied = copy_patent_svgs()
    print("copied svgs", copied)

    meta_path = ROOT / "metadata" / "variant_domain_metadata.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}

    # Prefer worker_1 (JAAROR lineage) — has positive fitness morphs
    w1 = ROOT / "outputs/rl_gym_data/worker_1/rl_training_dataset.jsonl"
    w0 = ROOT / "outputs/rl_gym_data/worker_0/rl_training_dataset.jsonl"
    rows1 = load_jsonl(w1) if w1.exists() else []
    rows0 = load_jsonl(w0) if w0.exists() else []

    best = None
    if rows1:
        out, best = fig_evolution(rows1, "JAAROR_lineage")
        print("evolution", out, "best", best["variant_id"], best["fitness"])
        coords = meta.get("NZ_JAAROR010000001.1_ORF_f1_65475_HIGH", {})
        print("hotspots", fig_mutation_map(rows1, coords, "JAAROR"))
    if rows0:
        out, _ = fig_evolution(rows0, "JAASWF_lineage")
        print("evolution", out)
        coords = meta.get("NZ_JAASWF010000012.1_ORF_f1_6252_HIGH", {})
        print("hotspots", fig_mutation_map(rows0, coords, "JAASWF"))

    # write wetlab shortlist table
    if rows1:
        top = sorted(rows1, key=lambda r: r["fitness"], reverse=True)[:8]
        lines = [
            "variant_id,baseline_id,generation,fitness,off_A,on_A,shift_A,iptm,n_mutations,mutations_head"
        ]
        for r in top:
            muts = r.get("mutations") or []
            lines.append(
                ",".join(
                    [
                        str(r.get("variant_id")),
                        str(r.get("baseline_id")),
                        str(r.get("generation")),
                        f"{r.get('fitness'):.4f}",
                        f"{r.get('off_dist_A'):.2f}",
                        f"{r.get('on_dist_A'):.2f}",
                        f"{(r.get('off_dist_A') or 0) - (r.get('on_dist_A') or 0):.2f}",
                        f"{r.get('iptm'):.4f}",
                        str(len(muts)),
                        ";".join(muts[:10]),
                    ]
                )
            )
        (FIG / "wetlab_shortlist_JAAROR.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
        print("wrote wetlab shortlist")


if __name__ == "__main__":
    main()
