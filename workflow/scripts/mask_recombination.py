import json
from pathlib import Path

from pathogenphyloflow.metrics import read_fasta


alignment_path = Path(snakemake.input.alignment)
gff_path = Path(snakemake.input.gff)
out_alignment = Path(snakemake.output.alignment)
out_summary = Path(snakemake.output.summary)

out_alignment.parent.mkdir(parents=True, exist_ok=True)

records = read_fasta(alignment_path)
lengths = {len(seq) for seq in records.values()}
if len(lengths) != 1:
    raise ValueError("Core alignment sequences do not all have the same length")
alignment_length = lengths.pop()

intervals = []
with open(gff_path, encoding="utf-8") as handle:
    for raw in handle:
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue
        try:
            start = max(1, int(fields[3]))
            end = min(alignment_length, int(fields[4]))
        except ValueError:
            continue
        if start <= end:
            intervals.append((start, end))

intervals.sort()
merged = []
for start, end in intervals:
    if not merged or start > merged[-1][1] + 1:
        merged.append([start, end])
    else:
        merged[-1][1] = max(merged[-1][1], end)

mask = [False] * alignment_length
for start, end in merged:
    for i in range(start - 1, end):
        mask[i] = True
masked_positions = sum(mask)

with open(out_alignment, "w", encoding="utf-8") as handle:
    for name, seq in records.items():
        chars = list(seq)
        for i, should_mask in enumerate(mask):
            if should_mask:
                chars[i] = "N"
        masked = "".join(chars)
        handle.write(f">{name}\n")
        for i in range(0, len(masked), 80):
            handle.write(masked[i : i + 80] + "\n")

summary = {
    "alignment_length": alignment_length,
    "raw_gff_intervals": len(intervals),
    "merged_intervals": len(merged),
    "masked_positions": masked_positions,
    "masked_fraction": 0.0 if alignment_length == 0 else masked_positions / alignment_length,
    "masking_strategy": "Global conservative mask of every alignment interval implicated by Gubbins",
}
out_summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
