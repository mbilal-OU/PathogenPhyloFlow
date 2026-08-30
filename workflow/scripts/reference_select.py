from __future__ import annotations

import csv
import json
import shutil
import statistics
import subprocess
import tempfile
from pathlib import Path

from workflow.lib.metrics import fasta_stats


assemblies = [Path(p) for p in snakemake.input.assemblies]
samples = list(snakemake.params.samples)
settings = dict(snakemake.params.settings)
out_ref = Path(snakemake.output.reference)
out_table = Path(snakemake.output.table)
out_json = Path(snakemake.output.summary)
log_path = Path(snakemake.log[0])

out_ref.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

if len(assemblies) != len(samples):
    raise ValueError("Reference selection received mismatched sample and assembly counts")

records = {}
for sample, assembly in zip(samples, assemblies):
    stats = fasta_stats(assembly)
    records[sample] = {
        "sample": sample,
        "assembly": str(assembly),
        **stats,
        "median_mash_distance": None,
        "eligible": True,
        "selected": False,
    }

mode = settings["mode"]
selection_reason = ""

if mode == "path":
    explicit = settings.get("path")
    if not explicit:
        raise ValueError("reference.mode=path requires reference.path")
    explicit_path = Path(explicit)
    if not explicit_path.exists():
        raise FileNotFoundError(f"Configured reference path does not exist: {explicit}")
    shutil.copyfile(explicit_path, out_ref)
    selected_sample = None
    selection_reason = "User supplied reference path"

elif mode == "sample":
    selected_sample = settings.get("sample")
    if selected_sample not in records:
        raise ValueError(f"Configured reference sample is not in the sample sheet: {selected_sample}")
    shutil.copyfile(records[selected_sample]["assembly"], out_ref)
    records[selected_sample]["selected"] = True
    selection_reason = "User selected reference sample"

elif mode == "auto":
    lengths = [r["total_length"] for r in records.values()]
    median_length = statistics.median(lengths)
    tolerance = float(settings["length_tolerance_fraction"])
    max_contigs = int(settings["max_contigs"])
    lower = median_length * (1 - tolerance)
    upper = median_length * (1 + tolerance)

    eligible = []
    for sample, record in records.items():
        record["eligible"] = (
            record["contigs"] <= max_contigs
            and lower <= record["total_length"] <= upper
        )
        if record["eligible"]:
            eligible.append(sample)
    if not eligible:
        eligible = list(samples)
        for sample in eligible:
            records[sample]["eligible"] = True

    with tempfile.TemporaryDirectory(prefix="pathogenphyloflow_mash_") as tempdir:
        prefix = Path(tempdir) / "candidates"
        sketch_cmd = ["mash", "sketch", "-o", str(prefix)] + [str(p) for p in assemblies]
        sketch = subprocess.run(sketch_cmd, capture_output=True, text=True)
        if sketch.returncode != 0:
            log_path.write_text(sketch.stdout + "\n" + sketch.stderr, encoding="utf-8")
            raise RuntimeError(f"Mash sketch failed; see {log_path}")
        dist = subprocess.run(
            ["mash", "dist", str(prefix) + ".msh", str(prefix) + ".msh"],
            capture_output=True,
            text=True,
        )
        log_path.write_text(
            sketch.stdout + "\n" + sketch.stderr + "\n" + dist.stdout + "\n" + dist.stderr,
            encoding="utf-8",
        )
        if dist.returncode != 0:
            raise RuntimeError(f"Mash distance calculation failed; see {log_path}")

    distances = {sample: [] for sample in samples}
    path_to_sample = {Path(p).stem: s for s, p in zip(samples, assemblies)}

    def sample_from_mash_name(value):
        stem = Path(value.split(":")[0]).stem
        return path_to_sample.get(stem, stem)

    for line in dist.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        a = sample_from_mash_name(fields[0])
        b = sample_from_mash_name(fields[1])
        if a in distances and b in distances and a != b:
            distances[a].append(float(fields[2]))

    for sample in samples:
        values = distances[sample]
        records[sample]["median_mash_distance"] = statistics.median(values) if values else 0.0

    selected_sample = min(
        eligible,
        key=lambda s: (
            records[s]["median_mash_distance"],
            -records[s]["n50"],
            records[s]["contigs"],
            s,
        ),
    )
    records[selected_sample]["selected"] = True
    shutil.copyfile(records[selected_sample]["assembly"], out_ref)
    selection_reason = (
        "Automatic medoid-style selection: eligible assembly with the lowest median Mash distance; "
        "N50 and contig count are deterministic tie breakers"
    )
else:
    raise ValueError(f"Unknown reference mode: {mode}")

fieldnames = [
    "sample",
    "assembly",
    "total_length",
    "contigs",
    "n50",
    "max_contig",
    "median_mash_distance",
    "eligible",
    "selected",
]
with open(out_table, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for sample in samples:
        writer.writerow(records[sample])

summary = {
    "mode": mode,
    "selected_sample": selected_sample,
    "selected_reference": str(out_ref),
    "reason": selection_reason,
    "candidates": records,
}
out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
