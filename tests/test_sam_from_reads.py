"""Tests for pipelines/shared/scripts/sam_from_reads.py — CIGAR parsing."""

import pytest
import sys
import os

# SamFromReads is a class with instance methods, import it
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipelines", "shared", "scripts"))
from sam_from_reads import SamFromReads


@pytest.fixture
def sam():
    """Create a SamFromReads instance for testing."""
    return SamFromReads()


class TestGetCigarLength:
    def test_simple_match(self, sam):
        assert sam.get_cigar_length("10M") == 10

    def test_multi_digit(self, sam):
        assert sam.get_cigar_length("150M") == 150

    def test_match_and_insertion(self, sam):
        # M and I both count toward query length
        assert sam.get_cigar_length("5M3I2M") == 10

    def test_deletion_excluded(self, sam):
        # D does not count toward query length
        assert sam.get_cigar_length("5M2D3M") == 8

    def test_complex_cigar(self, sam):
        # 10M + 5I + 3M = 18 (2D excluded)
        assert sam.get_cigar_length("10M2D5I3M") == 18

    def test_single_base(self, sam):
        assert sam.get_cigar_length("1M") == 1

    def test_soft_clip(self, sam):
        # S counts toward query length
        assert sam.get_cigar_length("5S10M5S") == 20
