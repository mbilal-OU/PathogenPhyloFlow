import json
import sys
from pathlib import Path


def _bootstrap_project_package():
    candidates = [Path.cwd(), *Path(__file__).resolve().parents]
    for root in candidates:
        if (root / "pathogenphyloflow" / "__init__.py").is_file():
            root_text = str(root)
            if root_text not in sys.path:
                sys.path.insert(0, root_text)
            return
    raise RuntimeError("Could not locate the PathogenPhyloFlow package from the Snakemake wrapper")


_bootstrap_project_package()
from pathogenphyloflow.masking import mask_alignment, parse_gubbins_intervals
from pathogenphyloflow.metrics import read_fasta


alignment_path = Path(snakemake.input.alignment)
gff_path = Path(snakemake.input.gff)
out_alignment = Path(snakemake.output.alignment)
out_summary = Path(snakemake.output.summary)
strategy = snakemake.config["recombination"].get("masking", "global")

out_alignment.parent.mkdir(parents=True, exist_ok=True)

records = read_fasta(alignment_path)
lengths = {len(seq) for seq in records.values()}
if len(lengths) != 1:
    raise ValueError("Core alignment sequences do not all have the same length")
alignment_length = lengths.pop()

intervals = parse_gubbins_intervals(gff_path, alignment_length)
masked_records, summary = mask_alignment(records, intervals, strategy=strategy)

with open(out_alignment, "w", encoding="utf-8") as handle:
    for name, seq in masked_records.items():
        handle.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            handle.write(seq[i : i + 80] + "\n")

out_summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
