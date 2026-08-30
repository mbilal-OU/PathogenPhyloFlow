from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

from workflow.lib.metrics import fasta_stats


row = dict(snakemake.params.row)
sample = row["sample"].strip()
assembly = row["assembly"].strip()
accession = row["accession"].strip()
out_fasta = Path(snakemake.output.assembly)
out_meta = Path(snakemake.output.metadata)
log_path = Path(snakemake.log[0])

out_fasta.parent.mkdir(parents=True, exist_ok=True)
out_meta.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

source = None
source_value = None

if assembly:
    source_path = Path(assembly)
    if not source_path.exists():
        raise FileNotFoundError(f"Assembly for {sample} does not exist: {assembly}")
    shutil.copyfile(source_path, out_fasta)
    source = "local_assembly"
    source_value = str(source_path)
else:
    if not (accession.startswith("GCA_") or accession.startswith("GCF_")):
        raise ValueError(
            f"Accession for {sample} must start with GCA_ or GCF_: {accession}"
        )
    if "REPLACE" in accession.upper():
        raise ValueError(f"Replace the placeholder accession for {sample}")
    with tempfile.TemporaryDirectory(prefix="pathogenphyloflow_") as tempdir:
        tempdir = Path(tempdir)
        archive = tempdir / "genome.zip"
        cmd = [
            "datasets",
            "download",
            "genome",
            "accession",
            accession,
            "--include",
            "genome",
            "--filename",
            str(archive),
        ]
        completed = subprocess.run(cmd, capture_output=True, text=True)
        log_path.write_text(completed.stdout + "\n" + completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(f"NCBI datasets failed for {sample}; see {log_path}")
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(tempdir / "download")
        fasta_candidates = list((tempdir / "download").rglob("*.fna"))
        if not fasta_candidates:
            raise RuntimeError(f"No genome FASTA found in NCBI download for {accession}")
        selected = max(fasta_candidates, key=lambda p: p.stat().st_size)
        shutil.copyfile(selected, out_fasta)
    source = "ncbi_accession"
    source_value = accession

stats = fasta_stats(out_fasta)
metadata = {
    "sample": sample,
    "source": source,
    "source_value": source_value,
    "date": row.get("date", ""),
    "location": row.get("location", ""),
    "host": row.get("host", ""),
    "assembly_stats": stats,
}
out_meta.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")

if assembly:
    log_path.write_text(f"Resolved local assembly {assembly} -> {out_fasta}\n", encoding="utf-8")
