"""Software-version provenance helpers.

Each workflow rule that shells out to an external bioinformatics tool records
that tool's version string into ``results/versions/<tool>.txt``. The
``collect_software_versions`` rule then aggregates those files into
``results/report/software_versions.json``. Version capture never fails a rule:
an unresolvable version is recorded as ``"unknown"``.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


def tool_version_string(tool: str, *args: str, timeout: int = 60) -> str:
    """Return the first line of ``tool *args`` output, or ``"unknown"``."""
    try:
        proc = subprocess.run(
            [tool, *args], capture_output=True, text=True, timeout=timeout
        )
        combined = (proc.stdout.strip() or proc.stderr.strip()).splitlines()
        return combined[0].strip() if combined else "unknown"
    except Exception:
        return "unknown"


def write_tool_version(path: str | Path, tool: str, *args: str) -> str:
    """Write ``tool``'s version string to ``path``; return the string."""
    version = tool_version_string(tool, *args)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(version + "\n", encoding="utf-8")
    return version
