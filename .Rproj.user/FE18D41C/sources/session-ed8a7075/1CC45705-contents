from kmer_analysis import validate_sequence
from kmer_analysis import update_kmer_count
from kmer_analysis import count_kmers_with_context


def test_validate_sequence_valid():
    assert validate_sequence("ATGC", 2) is True


def test_validate_sequence_too_short():
    assert validate_sequence("A", 2) is False


def test_validate_sequence_invalid_character():
    assert validate_sequence("ATGX", 2) is False


def test_update_kmer_count_new_kmer():
    result = update_kmer_count({}, "AT", "G")

    assert result["AT"]["count"] == 1
    assert result["AT"]["next_chars"]["G"] == 1


def test_count_kmers_with_context():
    result = count_kmers_with_context("ATGT", 2)

    assert result["AT"]["count"] == 1
    assert result["AT"]["next_chars"]["G"] == 1

    assert result["TG"]["count"] == 1
    assert result["TG"]["next_chars"]["T"] == 1

    assert result["GT"]["count"] == 1