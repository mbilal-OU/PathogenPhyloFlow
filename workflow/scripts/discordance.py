from __future__ import annotations

import csv
import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1]))
from workflow.lib.metrics import jaccard_distance


def read_snp_matrix(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("Empty SNP distance matrix")
    header = rows[0][1:]
    matrix = {}
    for row in rows[1:]:
        if not row:
            continue
        name = row[0]
        values = row[1:]
        matrix[name] = {other: int(float(value)) for other, value in zip(header, values)}
    return header, matrix


def read_accessory_sets(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        genes = {sample: set() for sample in samples}
        for row in reader:
            if not row:
                continue
            gene = row[0]
            for sample, value in zip(samples, row[1:]):
                try:
                    present = float(value) > 0
                except ValueError:
                    present = value.strip() not in {"", "0", "False", "false"}
                if present:
                    genes[sample].add(gene)
    return genes


snp_samples, snp = read_snp_matrix(snakemake.input.snp)
accessory = read_accessory_sets(snakemake.input.accessory)
common = [sample for sample in snp_samples if sample in accessory and sample in snp]
if len(common) < 2:
    raise ValueError("Fewer than two samples overlap between SNP and accessory matrices")

low_snp = int(snakemake.params.low_snp)
high_accessory = float(snakemake.params.high_accessory)
rows = []
for a, b in itertools.combinations(common, 2):
    snp_distance = snp[a][b]
    accessory_distance = jaccard_distance(accessory[a], accessory[b])
    flagged = snp_distance <= low_snp and accessory_distance >= high_accessory
    rows.append(
        {
            "sample_a": a,
            "sample_b": b,
            "snp_distance": snp_distance,
            "accessory_jaccard_distance": accessory_distance,
            "flag_low_snp_high_accessory": flagged,
        }
    )

rows.sort(key=lambda r: (not r["flag_low_snp_high_accessory"], r["snp_distance"], -r["accessory_jaccard_distance"]))
out_table = Path(snakemake.output.table)
out_table.parent.mkdir(parents=True, exist_ok=True)
with open(out_table, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

flagged = [row for row in rows if row["flag_low_snp_high_accessory"]]
summary = {
    "samples_compared": len(common),
    "pairs_compared": len(rows),
    "flagged_pairs": len(flagged),
    "low_snp_threshold": low_snp,
    "high_accessory_distance": high_accessory,
    "interpretation": "Flagged pairs are candidates for investigation, not transmission calls",
    "top_flagged_pairs": flagged[:20],
}
Path(snakemake.output.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
