"""Aggregate per-tool version files into results/report/software_versions.json.

Snakemake script (workflow/scripts/collect_versions.py). Inputs are the
``results/versions/<tool>.txt`` files written by each rule; outputs the JSON
catalogue plus interpreter/workflow metadata.
"""

import json
import platform
import subprocess
import sys
from pathlib import Path


def _git_commit() -> str | None:
    try:
        proc = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=15,
            cwd=str(Path(__file__).resolve().parents[2]),
        )
        commit = proc.stdout.strip()
        return commit if proc.returncode == 0 and commit else None
    except Exception:
        return None


version_files = [Path(p) for p in snakemake.input]
out_path = Path(snakemake.output.json)
out_path.parent.mkdir(parents=True, exist_ok=True)

tools = {}
for vf in version_files:
    text = vf.read_text(encoding="utf-8").strip() if vf.is_file() else "unknown"
    tools[vf.stem] = text.splitlines()[0] if text else "unknown"

catalogue = {
    "tools": dict(sorted(tools.items())),
    "python": platform.python_version(),
    "git_commit": _git_commit(),
}
try:
    import snakemake

    catalogue["snakemake"] = snakemake.__version__
except Exception:
    catalogue["snakemake"] = "unknown"

out_path.write_text(json.dumps(catalogue, indent=2, sort_keys=True) + "\n", encoding="utf-8")
