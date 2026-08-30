#!/usr/bin/env python3
"""Generate publication-style README figures from committed CI validation outputs."""

from __future__ import annotations

import csv
import json
from io import StringIO
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[2]
PREVIEW = ROOT / "examples" / "ecoli_o157" / "results_preview"
SAMPLES = ROOT / "examples" / "ecoli_o157" / "ci_samples.tsv"
OUT = ROOT / "docs" / "assets"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "savefig.bbox": "tight",
        "svg.fonttype": "none",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

summary = json.loads((PREVIEW / "validation_summary.json").read_text(encoding="utf-8"))
raw_newick = (PREVIEW / "raw.treefile").read_text(encoding="utf-8").strip()
filtered_newick = (PREVIEW / "recombination_filtered.treefile").read_text(encoding="utf-8").strip()

pairs = []
with (PREVIEW / "snp_accessory_discordance.tsv").open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        pairs.append(
            {
                "a": row["sample_a"],
                "b": row["sample_b"],
                "snp": int(row["snp_distance"]),
                "acc": float(row["accessory_jaccard_distance"]),
            }
        )

dates = {}
with SAMPLES.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row.get("date", "").strip():
            dates[row["sample"].strip()] = float(row["date"].strip())


def parse_tree(text: str, midpoint: bool = False):
    tree = Phylo.read(StringIO(text), "newick")
    if midpoint:
        tree.root_at_midpoint()
    return tree


def temporal_points():
    tree = parse_tree(filtered_newick, midpoint=True)
    values = []
    for tip in tree.get_terminals():
        if tip.name in dates:
            values.append((tip.name, dates[tip.name], tree.distance(tree.root, tip)))
    return sorted(values, key=lambda item: item[1])


points = temporal_points()
x_year = np.array([item[1] for item in points], dtype=float)
y_rtt = np.array([item[2] for item in points], dtype=float)
slope, intercept = np.polyfit(x_year, y_rtt, 1)
predicted = slope * x_year + intercept
r2 = 1 - np.sum((y_rtt - predicted) ** 2) / np.sum((y_rtt - np.mean(y_rtt)) ** 2)

# Figure 1: four-panel validation overview.
fig = plt.figure(figsize=(10.5, 7.5), constrained_layout=True)
gs = fig.add_gridspec(2, 2)

ax = fig.add_subplot(gs[0, 0])
masked_fraction = float(summary["recombination"]["masked_fraction"])
ax.barh([0], [1 - masked_fraction], height=0.35, label="Retained")
ax.barh([0], [masked_fraction], left=[1 - masked_fraction], height=0.35, label="Masked")
ax.set_xlim(0, 1)
ax.set_ylim(-0.7, 0.7)
ax.set_yticks([])
ax.set_xlabel("Fraction of core alignment")
ax.set_title("A  Recombination masking", loc="left", fontweight="bold")
ax.legend(frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0))
ax.text(
    0.0,
    -0.38,
    f"Alignment length = {summary['recombination']['alignment_length']:,} bp\n"
    f"Masked = {summary['recombination']['masked_positions']:,} bp ({masked_fraction * 100:.2f}%)\n"
    f"Merged recombinant intervals = {summary['recombination']['merged_intervals']}",
    fontsize=8,
    va="top",
)
ax.spines["left"].set_visible(False)

ax = fig.add_subplot(gs[0, 1])
raw_length = float(summary["tree_comparison"]["raw_total_branch_length"])
filtered_length = float(summary["tree_comparison"]["filtered_total_branch_length"])
ratio = float(summary["tree_comparison"]["branch_length_ratio_filtered_to_raw"])
ax.bar(["Raw", "Filtered"], [raw_length, filtered_length], width=0.55)
ax.set_ylabel("Total branch length")
ax.ticklabel_format(axis="y", style="sci", scilimits=(-4, -4))
ax.set_title("B  Tree comparison", loc="left", fontweight="bold")
ax.text(
    0.98,
    0.96,
    f"Normalized RF = {summary['tree_comparison']['normalized_robinson_foulds']:.0f}\nFiltered/raw = {ratio:.3f}",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=8,
)
ax.grid(axis="y", alpha=0.12)

ax = fig.add_subplot(gs[1, 0])
for row in pairs:
    ax.scatter(row["snp"], row["acc"], s=28)
low_snp = float(summary["snp_accessory_discordance"]["low_snp_threshold"])
high_acc = float(summary["snp_accessory_discordance"]["high_accessory_distance"])
ax.axvline(low_snp, linestyle="--", linewidth=1)
ax.axhline(high_acc, linestyle="--", linewidth=1)
ax.set_xlim(0, 740)
ax.set_ylim(0, 0.30)
ax.set_xlabel("Core-SNP distance")
ax.set_ylabel("Accessory Jaccard distance")
ax.set_title("C  SNP/accessory comparison", loc="left", fontweight="bold")
ax.text(
    0.98,
    0.05,
    f"{summary['snp_accessory_discordance']['pairs_compared']} pairwise comparisons\n"
    f"{summary['snp_accessory_discordance']['flagged_pairs']} flagged pairs\n"
    f"{summary['accessory_gene_families']:,} accessory gene families",
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    fontsize=8,
)
ax.grid(alpha=0.12)

ax = fig.add_subplot(gs[1, 1])
ax.scatter(x_year, y_rtt, s=30)
for name, year, distance in points:
    ax.annotate(name, (year, distance), xytext=(4, 4), textcoords="offset points", fontsize=7)
xx = np.linspace(x_year.min() - 1, x_year.max() + 1, 100)
ax.plot(xx, slope * xx + intercept, linestyle="--", linewidth=1.2)
ax.set_xlabel("Sampling year")
ax.set_ylabel("Root-to-tip distance")
ax.ticklabel_format(axis="y", style="sci", scilimits=(-5, -5))
ax.set_title("D  Temporal diagnostic", loc="left", fontweight="bold")
ax.text(
    0.02,
    0.04,
    f"R² = {r2:.4f}; slope < 0\nScreen not passed",
    transform=ax.transAxes,
    fontsize=8,
    va="bottom",
)
ax.grid(alpha=0.12)

fig.suptitle(
    "PathogenPhyloFlow CI validation: E. coli O157:H7 four-genome subset",
    fontsize=12,
    fontweight="bold",
)
fig.savefig(OUT / "ecoli_o157_validation_snapshot.svg")
plt.close(fig)

# Figure 2: actual raw and recombination-filtered trees.
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), constrained_layout=True)
for ax, tree, title in zip(
    axes,
    [parse_tree(raw_newick), parse_tree(filtered_newick)],
    ["Raw core-SNP tree", "Recombination-filtered tree"],
):
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=True)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("Substitutions per site")
    ax.set_ylabel("")
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.ticklabel_format(axis="x", style="sci", scilimits=(-4, -4))
    ax.grid(axis="x", alpha=0.12)
fig.suptitle("Recombination-aware phylogeny comparison", fontsize=12, fontweight="bold")
fig.savefig(OUT / "ecoli_o157_tree_comparison.svg")
plt.close(fig)

# Figure 3: SNP versus accessory-genome distance with exact configured thresholds.
fig, ax = plt.subplots(figsize=(7.4, 5.0), constrained_layout=True)
label_offsets = {
    ("EC4115", "TW14359"): (7, 4),
    ("EDL933", "Sakai"): (7, 6),
    ("EC4115", "Sakai"): (-85, 7),
    ("Sakai", "TW14359"): (7, -15),
    ("EC4115", "EDL933"): (-95, -13),
    ("EDL933", "TW14359"): (7, 6),
}
for row in pairs:
    ax.scatter(row["snp"], row["acc"], s=34, zorder=3)
    dx, dy = label_offsets[(row["a"], row["b"])]
    ax.annotate(
        f"{row['a']}–{row['b']}",
        (row["snp"], row["acc"]),
        xytext=(dx, dy),
        textcoords="offset points",
        fontsize=7,
    )
ax.axvline(low_snp, linestyle="--", linewidth=1)
ax.axhline(high_acc, linestyle="--", linewidth=1)
ax.set_xlim(0, 740)
ax.set_ylim(0, 0.30)
ax.set_xlabel("Pairwise core-SNP distance")
ax.set_ylabel("Accessory-genome Jaccard distance")
ax.set_title("Core-SNP versus accessory-genome distance", loc="left", fontweight="bold")
ax.text(low_snp + 8, 0.286, "Low-SNP threshold = 10", fontsize=8, va="top")
ax.text(735, high_acc + 0.005, "High accessory-distance threshold = 0.25", fontsize=8, ha="right", va="bottom")
ax.grid(alpha=0.12)
fig.savefig(OUT / "ecoli_o157_snp_accessory.svg")
plt.close(fig)

# Figure 4: midpoint-rooted temporal diagnostic using the same calculation as the workflow.
fig, ax = plt.subplots(figsize=(7.4, 5.0), constrained_layout=True)
ax.scatter(x_year, y_rtt, s=38, zorder=3)
for name, year, distance in points:
    ax.annotate(name, (year, distance), xytext=(7, 5), textcoords="offset points", fontsize=8)
ax.plot(xx, slope * xx + intercept, linestyle="--", linewidth=1.4)
ax.set_xlabel("Sampling year")
ax.set_ylabel("Root-to-tip genetic distance")
ax.set_title("Root-to-tip temporal diagnostic", loc="left", fontweight="bold")
ax.ticklabel_format(axis="y", style="sci", scilimits=(-5, -5))
ax.grid(alpha=0.12)
ax.text(
    0.02,
    0.04,
    f"n = {len(points)}; span = {summary['temporal']['sampling_span_years']:.0f} years; "
    f"slope = {slope:.2e}; R² = {r2:.4f}\nTemporal screen not passed",
    transform=ax.transAxes,
    fontsize=8,
    va="bottom",
)
fig.savefig(OUT / "ecoli_o157_temporal_screen.svg")
plt.close(fig)

print("Generated academic validation figures in docs/assets")
