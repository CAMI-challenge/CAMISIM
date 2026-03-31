"""Tests for pipelines/shared/scripts/sam_from_reads.py."""

import io
import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipelines", "shared", "scripts"))
from sam_from_reads import SamFromReads


@pytest.fixture
def sam():
    """Create a SamFromReads instance for testing."""
    return SamFromReads()


# ---------------------------------------------------------------------------
# get_cigar_length
# ---------------------------------------------------------------------------
class TestGetCigarLength:
    def test_simple_match(self, sam):
        assert sam.get_cigar_length("10M") == 10

    def test_multi_digit(self, sam):
        assert sam.get_cigar_length("150M") == 150

    def test_match_and_insertion(self, sam):
        assert sam.get_cigar_length("5M3I2M") == 10

    def test_deletion_excluded(self, sam):
        assert sam.get_cigar_length("5M2D3M") == 8

    def test_complex_cigar(self, sam):
        assert sam.get_cigar_length("10M2D5I3M") == 18

    def test_single_base(self, sam):
        assert sam.get_cigar_length("1M") == 1

    def test_soft_clip(self, sam):
        assert sam.get_cigar_length("5S10M5S") == 20

    def test_only_insertion(self, sam):
        assert sam.get_cigar_length("5I") == 5

    def test_only_deletion(self, sam):
        assert sam.get_cigar_length("5D") == 0

    def test_large_numbers(self, sam):
        assert sam.get_cigar_length("1000M50I20D30M") == 1080


# ---------------------------------------------------------------------------
# get_cigars_nanosim
# ---------------------------------------------------------------------------
class TestGetCigarsNanosim:
    def _write_error_profile(self, tmp_path, lines):
        ep = tmp_path / "error_profile.txt"
        ep.write_text("\n".join(lines) + "\n")
        return str(ep)

    def test_insertion_only(self, sam, tmp_path):
        lines = [
            "Seq\tPos\tType\tLen\tRef\tQuery",
            "refA_100_aligned_0_F_0_500_0\t10\tins\t3\tACG\tACGTTT",
        ]
        ep = self._write_error_profile(tmp_path, lines)
        cigars = sam.get_cigars_nanosim(ep)
        assert "refA-100-0" in cigars
        cigar_str, _ = cigars["refA-100-0"]
        assert "10M" in cigar_str
        assert "3I" in cigar_str

    def test_deletion_only(self, sam, tmp_path):
        lines = [
            "Seq\tPos\tType\tLen\tRef\tQuery",
            "refA_100_aligned_0_F_0_500_0\t10\tdel\t2\tAC\t--",
        ]
        ep = self._write_error_profile(tmp_path, lines)
        cigars = sam.get_cigars_nanosim(ep)
        cigar_str, _ = cigars["refA-100-0"]
        assert "10M" in cigar_str
        assert "2D" in cigar_str

    def test_mismatches_ignored(self, sam, tmp_path):
        lines = [
            "Seq\tPos\tType\tLen\tRef\tQuery",
            "refA_100_aligned_0_F_0_500_0\t10\tmis\t1\tA\tG",
        ]
        ep = self._write_error_profile(tmp_path, lines)
        cigars = sam.get_cigars_nanosim(ep)
        assert len(cigars) == 0  # mismatches don't produce entries

    def test_multiple_errors_same_read(self, sam, tmp_path):
        lines = [
            "Seq\tPos\tType\tLen\tRef\tQuery",
            "refA_100_aligned_0_F_0_500_0\t10\tins\t2\tAC\tACGG",
            "refA_100_aligned_0_F_0_500_0\t20\tdel\t3\tACG\t---",
        ]
        ep = self._write_error_profile(tmp_path, lines)
        cigars = sam.get_cigars_nanosim(ep)
        cigar_str, _ = cigars["refA-100-0"]
        assert "I" in cigar_str
        assert "D" in cigar_str
