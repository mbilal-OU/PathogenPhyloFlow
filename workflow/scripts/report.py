from __future__ import annotations

import csv
import html
import json
import subprocess
from pathlib import Path


def load_json(path, default=None):
    path = Path(path)
    if not path.exists():
        return default if default is not None else {}
    return json.loads(path.read_text(encoding="utf-8"))


def count_tsv_data_rows(path):
    path = Path(path)
    if not path.exists():
        return 0
    with open(path, encoding="utf-8") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return None


with open(snakemake.params.samples, newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    sample_rows = list(reader)
    sample_fields = list(reader.fieldnames or [])

reference = load_json("results/reference/reference_selection.json")
mask = load_json("results/recombination/mask_summary.json", {})
tree_compare = load_json("results/phylogeny/tree_comparison.json", {})
temporal = load_json("results/temporal/temporal_screen.json", {})
discordance = load_json("results/integration/snp_accessory_summary.json", {})

accessory_gene_count = count_tsv_data_rows("results/accessory/panaroo/gene_presence_absence.Rtab")
resfinder_rows = count_tsv_data_rows("results/functional/resfinder_summary.tsv")
vfdb_rows = count_tsv_data_rows("results/functional/vfdb_summary.tsv")

report_metadata_fields = [
    field for field in sample_fields
    if field not in {"assembly", "accession"}
]

summary = {
    "project": "PathogenPhyloFlow",
    "samples": len(sample_rows),
    "sample_metadata_fields": report_metadata_fields,
    "commit": git_commit(),
    "reference": reference,
    "recombination": mask,
    "tree_comparison": tree_compare,
    "accessory_gene_families": accessory_gene_count,
    "snp_accessory_discordance": discordance,
    "temporal": temporal,
    "functional": {
        "enabled": bool(snakemake.params.functional),
        "resfinder_summary_rows": resfinder_rows,
        "vfdb_summary_rows": vfdb_rows,
    },
    "guardrails": [
        "SNP proximity is not treated as proof of direct transmission.",
        "Recombination effects are retained and reported when enabled.",
        "Root-to-tip regression is treated as a temporal diagnostic, not proof of a molecular clock.",
        "Accessory-genome discordance flags are investigation targets, not conclusions.",
    ],
}

out_html = Path(snakemake.output.html)
out_json = Path(snakemake.output.summary)
out_config = Path(snakemake.output.config)
out_html.parent.mkdir(parents=True, exist_ok=True)
out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
out_config.write_text(json.dumps(dict(snakemake.config), indent=2, sort_keys=True) + "\n", encoding="utf-8")

selected = reference.get("selected_sample") or "user-supplied reference"
masked_fraction = mask.get("masked_fraction")
masked_text = "not run" if masked_fraction is None else f"{masked_fraction:.2%}"
rf = tree_compare.get("normalized_robinson_foulds")
rf_text = "not available" if rf is None else f"{rf:.3f}"
temporal_status = temporal.get("treetime_status", "not available")
r2 = temporal.get("r2")
r2_text = "not available" if r2 is None else f"{r2:.3f}"
flagged = discordance.get("flagged_pairs", 0)

cards = [
    ("Samples", len(sample_rows), "Resolved pathogen genomes"),
    ("Reference", selected, "Recorded selection rationale"),
    ("Recombination masked", masked_text, "Fraction of core alignment"),
    ("Tree change", rf_text, "Normalized RF distance"),
    ("Accessory families", accessory_gene_count if accessory_gene_count else "off", "Panaroo gene families"),
    ("Discordant pairs", flagged, "Low SNP, high accessory distance"),
    ("Temporal screen R²", r2_text, temporal_status),
]

card_html = "".join(
    f"<div class='card'><div class='label'>{html.escape(str(label))}</div>"
    f"<div class='value'>{html.escape(str(value))}</div>"
    f"<div class='note'>{html.escape(str(note))}</div></div>"
    for label, value, note in cards
)

guardrails = "".join(f"<li>{html.escape(item)}</li>" for item in summary["guardrails"])

header_html = "".join(
    f"<th>{html.escape(field.replace('_', ' ').title())}</th>"
    for field in report_metadata_fields
)
rows_html = "".join(
    "<tr>"
    + "".join(
        f"<td>{html.escape(row.get(field, '') or '')}</td>"
        for field in report_metadata_fields
    )
    + "</tr>"
    for row in sample_rows[:200]
)

page = f"""<!doctype html>
<html lang='en'>
<head>
<meta charset='utf-8'>
<meta name='viewport' content='width=device-width, initial-scale=1'>
<title>PathogenPhyloFlow report</title>
<style>
:root {{ color-scheme: light; --ink:#17212b; --muted:#66717d; --line:#d9e0e6; --panel:#f7f9fb; --accent:#285a70; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; font-family:Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; color:var(--ink); background:white; }}
main {{ max-width:1180px; margin:0 auto; padding:44px 24px 80px; }}
h1 {{ margin:0; font-size:36px; letter-spacing:-0.03em; }}
.subtitle {{ margin-top:8px; color:var(--muted); max-width:820px; line-height:1.55; }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:12px; margin:30px 0; }}
.card {{ border:1px solid var(--line); border-radius:12px; padding:18px; background:var(--panel); }}
.label {{ font-size:12px; text-transform:uppercase; letter-spacing:.08em; color:var(--muted); }}
.value {{ font-size:27px; font-weight:700; margin-top:6px; color:var(--accent); overflow-wrap:anywhere; }}
.note {{ font-size:12px; color:var(--muted); margin-top:5px; }}
section {{ margin-top:34px; }}
h2 {{ font-size:20px; margin-bottom:12px; }}
.notice {{ border-left:4px solid var(--accent); background:var(--panel); padding:14px 18px; line-height:1.55; }}
.table-wrap {{ width:100%; overflow-x:auto; }}
table {{ width:100%; border-collapse:collapse; font-size:14px; white-space:nowrap; }}
th, td {{ text-align:left; border-bottom:1px solid var(--line); padding:10px 8px; vertical-align:top; }}
th {{ color:var(--muted); font-weight:600; }}
code {{ background:var(--panel); padding:2px 5px; border-radius:4px; }}
.footer {{ color:var(--muted); font-size:12px; margin-top:46px; }}
</style>
</head>
<body><main>
<h1>PathogenPhyloFlow</h1>
<p class='subtitle'>Integrated pathogen genomic epidemiology report. Core SNP evolution, recombination, accessory-genome change, functional context, and temporal diagnostics are kept as separate evidence layers before interpretation.</p>
<div class='grid'>{card_html}</div>
<section><h2>Interpretation guardrails</h2><div class='notice'><ul>{guardrails}</ul></div></section>
<section><h2>Samples</h2><div class='table-wrap'><table><thead><tr>{header_html}</tr></thead><tbody>{rows_html}</tbody></table></div></section>
<section><h2>Stable outputs</h2><p><code>results/phylogeny/final.treefile</code> &nbsp; <code>results/phylogeny/final.snp_dist.tsv</code> &nbsp; <code>results/integration/snp_accessory_discordance.tsv</code></p></section>
<p class='footer'>Generated by PathogenPhyloFlow. Commit: {html.escape(summary['commit'] or 'unavailable')}</p>
</main></body></html>
"""
out_html.write_text(page, encoding="utf-8")
