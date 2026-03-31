"""Tests for pipelines/metagenomic/scripts/wgsim_to_sam.py — CIGAR generation."""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipelines", "metagenomic", "scripts"))
from wgsim_to_sam import get_cigar


class TestGetCigar:
    def test_all_matches(self):
        assert get_cigar("ACGT", "ACGT") == "4M"

    def test_mismatch_counted_as_match(self):
        # mismatches are represented by M in CIGAR
        assert get_cigar("ACGT", "AGGT") == "4M"

    def test_deletion_in_read(self):
        assert get_cigar("ACGT", "A-GT") == "1M1D2M"

    def test_insertion_in_read(self):
        assert get_cigar("A-GT", "ACGT") == "1M1I2M"

    def test_consecutive_deletions(self):
        assert get_cigar("ACGT", "A--T") == "1M2D1M"

    def test_consecutive_insertions(self):
        assert get_cigar("A--T", "ACGT") == "1M2I1M"

    def test_single_base(self):
        assert get_cigar("A", "A") == "1M"

    def test_complex_alignment(self):
        # 2M + 1D + 1I + 2M (zip pairs: AC match, G/- is D, -/A is I, TT match)
        assert get_cigar("ACG-TT", "AC-ATT") == "2M1D1I2M"

    def test_all_deletions(self):
        assert get_cigar("AAAA", "----") == "4D"

    def test_all_insertions(self):
        assert get_cigar("----", "AAAA") == "4I"
