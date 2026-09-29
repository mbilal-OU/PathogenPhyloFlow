"""Unit tests for pathogenphyloflow.masking (recombination masking strategies)."""

import pytest

from pathogenphyloflow.masking import (
    STRATEGIES,
    mask_alignment,
    merge_intervals,
    parse_gubbins_intervals,
)


@pytest.fixture()
def records():
    return {
        "taxonA": "ACGTACGTACGTACGTACGT",
        "taxonB": "ACGTACGTACGTACGTACGT",
        "taxonC": "ACGTACGTACGTACGTACGT",
    }


@pytest.fixture()
def intervals():
    # (seqname, start, end), 1-based inclusive
    return [
        ("taxonA", 5, 8),    # 4 positions in taxonA
        ("taxonB", 10, 12),  # 3 positions in taxonB
        ("ghost", 1, 3),     # unknown taxon: skipped under per_taxon
    ]


def test_merge_intervals_overlapping_and_adjacent():
    assert merge_intervals([(5, 8), (7, 10), (12, 14), (15, 15)]) == [
        (5, 10),
        (12, 15),
    ]


def test_parse_gubbins_intervals_skips_comments_and_clamps(tmp_path):
    gff = tmp_path / "rec.gff"
    gff.write_text(
        "##gff-version 3\n"
        "taxonA\tgubbins\trecombination\t5\t8\t.\t.\t.\tID=1\n"
        "taxonB\tgubbins\trecombination\t0\t25\t.\t.\t.\tID=2\n"  # clamped to 1..20
        "badline\n",
        encoding="utf-8",
    )
    parsed = parse_gubbins_intervals(gff, 20)
    assert parsed == [("taxonA", 5, 8), ("taxonB", 1, 20)]


def test_global_strategy_masks_union_in_every_sequence(records, intervals):
    masked, summary = mask_alignment(records, intervals, strategy="global")
    # union of intervals = 1-3, 5-8, 10-12 -> 10 positions in every sequence
    for name, seq in masked.items():
        assert seq[0:3] == "NNN"      # ghost interval still masks globally
        assert seq[4:8] == "NNNN"
        assert seq[9:12] == "NNN"
        assert seq[3] == "T" and seq[8] == "A" and seq[12] == "A"
    assert summary["masking_strategy"] == "global"
    assert summary["masked_positions"] == 30
    assert summary["masked_positions_per_taxon"] == {
        "taxonA": 10,
        "taxonB": 10,
        "taxonC": 10,
    }
    assert summary["merged_intervals"] == 3
    assert summary["masked_fraction"] == pytest.approx(30 / 60)


def test_per_taxon_strategy_masks_only_implicated_taxa(records, intervals):
    masked, summary = mask_alignment(records, intervals, strategy="per_taxon")
    assert masked["taxonA"][4:8] == "NNNN"
    assert masked["taxonA"][9:12] != "NNN"  # taxonB's interval not applied
    assert masked["taxonB"][9:12] == "NNN"
    assert masked["taxonB"][4:8] != "NNNN"
    assert set(masked["taxonC"]) == {"A", "C", "G", "T"}  # untouched
    assert summary["masking_strategy"] == "per_taxon"
    assert summary["masked_positions"] == 7
    assert summary["masked_positions_per_taxon"] == {
        "taxonA": 4,
        "taxonB": 3,
        "taxonC": 0,
    }
    assert summary["skipped_unknown_taxa_intervals"] == 1


def test_unknown_strategy_raises(records, intervals):
    with pytest.raises(ValueError, match="Unknown masking strategy"):
        mask_alignment(records, intervals, strategy="sideways")


def test_strategies_documented():
    assert set(STRATEGIES) == {"global", "per_taxon"}
