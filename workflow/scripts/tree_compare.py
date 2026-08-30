from __future__ import annotations

import csv
import json
from pathlib import Path

from Bio import Phylo


def split_set(tree):
    tips = {tip.name for tip in tree.get_terminals()}
    splits = set()
    for clade in tree.get_nonterminals(order="postorder"):
        subset = frozenset(t.name for t in clade.get_terminals())
        other = frozenset(tips - subset)
        if len(subset) <= 1 or len(other) <= 1:
            continue
        canonical = subset if (len(subset), sorted(subset)) <= (len(other), sorted(other)) else other
        splits.add(canonical)
    return tips, splits


def total_branch_length(tree):
    return sum((clade.branch_length or 0.0) for clade in tree.find_clades())


raw = Phylo.read(snakemake.input.raw, "newick")
filtered = Phylo.read(snakemake.input.filtered, "newick")
raw_tips, raw_splits = split_set(raw)
filtered_tips, filtered_splits = split_set(filtered)
if raw_tips != filtered_tips:
    missing_raw = sorted(filtered_tips - raw_tips)
    missing_filtered = sorted(raw_tips - filtered_tips)
    raise ValueError(
        f"Tree tip sets differ. Missing from raw: {missing_raw}; missing from filtered: {missing_filtered}"
    )

rf = len(raw_splits.symmetric_difference(filtered_splits))
denominator = len(raw_splits) + len(filtered_splits)
normalized_rf = 0.0 if denominator == 0 else rf / denominator
raw_length = total_branch_length(raw)
filtered_length = total_branch_length(filtered)

summary = {
    "tips": len(raw_tips),
    "raw_internal_splits": len(raw_splits),
    "filtered_internal_splits": len(filtered_splits),
    "robinson_foulds_symmetric_difference": rf,
    "normalized_robinson_foulds": normalized_rf,
    "raw_total_branch_length": raw_length,
    "filtered_total_branch_length": filtered_length,
    "branch_length_ratio_filtered_to_raw": None if raw_length == 0 else filtered_length / raw_length,
}

Path(snakemake.output.json).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.output.json).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
with open(snakemake.output.tsv, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["metric", "value"])
    for key, value in summary.items():
        writer.writerow([key, value])
