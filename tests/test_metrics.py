from pathlib import Path

import pytest

from pathogenphyloflow.metrics import decimal_year, fasta_stats, jaccard_distance, linear_regression, n50


def test_n50():
    # Total length is 20; the first 10-bp contig reaches 50% exactly.
    assert n50([10, 9, 1]) == 10
    # Here the first 9-bp contig is below 50%; the second 9-bp contig crosses it.
    assert n50([9, 9, 2]) == 9
    assert n50([]) == 0


def test_jaccard_distance():
    assert jaccard_distance({"a", "b"}, {"a", "b"}) == 0.0
    assert jaccard_distance({"a"}, {"b"}) == 1.0
    assert jaccard_distance(set(), set()) == 0.0


def test_decimal_year_ordering():
    assert decimal_year("2025-01-01") < decimal_year("2025-12-31")
    assert 2025.0 <= decimal_year("2025-06") < 2026.0


def test_linear_regression_perfect_fit():
    result = linear_regression([1.0, 2.0, 3.0], [2.0, 4.0, 6.0])
    assert result["slope"] == pytest.approx(2.0)
    assert result["r2"] == pytest.approx(1.0)


def test_fasta_stats(tmp_path: Path):
    fasta = tmp_path / "x.fna"
    fasta.write_text(">a\nAAAAAA\n>b\nAAAA\n", encoding="utf-8")
    stats = fasta_stats(fasta)
    assert stats == {"contigs": 2, "total_length": 10, "n50": 6, "max_contig": 6}
