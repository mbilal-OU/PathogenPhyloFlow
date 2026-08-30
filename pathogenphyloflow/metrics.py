from __future__ import annotations

from datetime import date
from pathlib import Path
from typing import Iterable


def read_fasta(path: str | Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current:
                    raise ValueError(f"Empty FASTA identifier in {path}")
                if current in records:
                    raise ValueError(f"Duplicate FASTA identifier {current!r} in {path}")
                records[current] = []
            elif current is None:
                raise ValueError(f"Sequence found before FASTA header in {path}")
            else:
                records[current].append(line)
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return {name: "".join(parts) for name, parts in records.items()}


def n50(lengths: Iterable[int]) -> int:
    values = sorted((int(x) for x in lengths if int(x) > 0), reverse=True)
    if not values:
        return 0
    half = sum(values) / 2
    cumulative = 0
    for value in values:
        cumulative += value
        if cumulative >= half:
            return value
    return 0


def fasta_stats(path: str | Path) -> dict[str, int]:
    records = read_fasta(path)
    lengths = [len(seq) for seq in records.values()]
    return {
        "contigs": len(lengths),
        "total_length": sum(lengths),
        "n50": n50(lengths),
        "max_contig": max(lengths),
    }


def jaccard_distance(a: set[str], b: set[str]) -> float:
    union = a | b
    if not union:
        return 0.0
    return 1.0 - (len(a & b) / len(union))


def decimal_year(value: str) -> float:
    value = value.strip()
    if not value:
        raise ValueError("Empty date")
    parts = value.split("-")
    year = int(parts[0])
    if len(parts) == 1:
        return float(year) + 0.5
    month = int(parts[1])
    day = int(parts[2]) if len(parts) >= 3 else 15
    current = date(year, month, day)
    start = date(year, 1, 1)
    end = date(year + 1, 1, 1)
    return year + ((current - start).days / (end - start).days)


def linear_regression(x: list[float], y: list[float]) -> dict[str, float]:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("Regression requires equal-length x and y with at least two values")
    xbar = sum(x) / len(x)
    ybar = sum(y) / len(y)
    ssx = sum((v - xbar) ** 2 for v in x)
    ssy = sum((v - ybar) ** 2 for v in y)
    if ssx == 0:
        return {"slope": 0.0, "intercept": ybar, "r2": 0.0}
    cross = sum((a - xbar) * (b - ybar) for a, b in zip(x, y))
    slope = cross / ssx
    intercept = ybar - slope * xbar
    r2 = 0.0 if ssy == 0 else (cross * cross) / (ssx * ssy)
    return {"slope": slope, "intercept": intercept, "r2": max(0.0, min(1.0, r2))}
