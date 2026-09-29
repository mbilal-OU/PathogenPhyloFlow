"""Recombination masking strategies for core-genome alignments.

Two strategies are supported, selected with ``recombination.masking``:

- ``"global"`` (default): every alignment interval implicated by Gubbins (the
  union over all taxa) is masked with ``N`` in *every* sequence. This is
  conservative — it can only remove phylogenetic signal, never create false
  signal — but it discards informative sites for taxa that were not involved
  in a recombination event.
- ``"per_taxon"``: each Gubbins interval is masked only in the taxon named in
  the GFF ``seqname`` column, preserving signal for uninvolved taxa. This
  matches what Gubbins itself does when it builds its own filtered alignment
  (``gubbins.filtered_polymorphic_sites.fasta``).

The default is ``"global"`` to preserve the behavior the workflow was
validated against; ``"per_taxon"`` is available for analyses where retaining
per-taxon signal matters more than a conservative mask.
"""

from __future__ import annotations

from pathlib import Path

STRATEGIES = ("global", "per_taxon")


def parse_gubbins_intervals(gff_path: str | Path, alignment_length: int):
    """Parse a Gubbins ``recombination_predictions.gff`` file.

    Returns a list of ``(seqname, start, end)`` tuples with 1-based inclusive
    coordinates clamped to ``[1, alignment_length]``. Comment and malformed
    lines are skipped.
    """
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
                intervals.append((fields[0], start, end))
    return intervals


def merge_intervals(intervals):
    """Merge a list of ``(start, end)`` 1-based inclusive intervals."""
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(s, e) for s, e in merged]


def mask_alignment(records: dict[str, str], intervals, strategy: str = "global"):
    """Mask recombinant intervals with ``N``.

    ``records`` maps sequence name to sequence (all must share one length);
    ``intervals`` is a list of ``(seqname, start, end)`` 1-based inclusive
    tuples as returned by :func:`parse_gubbins_intervals`.

    Returns ``(masked_records, summary)``. Intervals naming a sequence absent
    from ``records`` are skipped and counted in
    ``summary["skipped_unknown_taxa_intervals"]`` (per-taxon strategy only).
    """
    if strategy not in STRATEGIES:
        raise ValueError(
            f"Unknown masking strategy {strategy!r}; expected one of {STRATEGIES}"
        )
    lengths = {len(seq) for seq in records.values()}
    if len(lengths) != 1:
        raise ValueError("Alignment sequences do not all have the same length")
    alignment_length = lengths.pop() if lengths else 0

    merged_count = 0
    skipped_unknown = 0
    if strategy == "global":
        union = merge_intervals([(s, e) for _, s, e in intervals])
        merged_count = len(union)
        per_taxon = {name: union for name in records}
    else:
        by_taxon: dict[str, list] = {}
        for seqname, start, end in intervals:
            if seqname not in records:
                skipped_unknown += 1
                continue
            by_taxon.setdefault(seqname, []).append((start, end))
        per_taxon = {}
        for name in records:
            merged_taxon = merge_intervals(by_taxon.get(name, []))
            merged_count += len(merged_taxon)
            per_taxon[name] = merged_taxon

    masked_records = {}
    masked_per_taxon = {}
    for name, seq in records.items():
        chars = list(seq)
        count = 0
        for start, end in per_taxon[name]:
            for i in range(start - 1, min(end, alignment_length)):
                if chars[i] != "N":
                    count += 1
                chars[i] = "N"
        masked_records[name] = "".join(chars)
        masked_per_taxon[name] = count

    total_masked = sum(masked_per_taxon.values())
    total_cells = alignment_length * len(records)
    summary = {
        "alignment_length": alignment_length,
        "masking_strategy": strategy,
        "strategy_description": (
            "Global conservative mask of every alignment interval implicated by Gubbins"
            if strategy == "global"
            else "Per-taxon mask: each Gubbins interval masked only in its implicated taxon"
        ),
        "raw_gff_intervals": len(intervals),
        "merged_intervals": merged_count,
        "masked_positions": total_masked,
        "masked_fraction": 0.0 if total_cells == 0 else total_masked / total_cells,
        "masked_positions_per_taxon": masked_per_taxon,
    }
    if strategy == "per_taxon":
        summary["skipped_unknown_taxa_intervals"] = skipped_unknown
    return masked_records, summary
