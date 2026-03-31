"""Tests for pipelines/metatranscriptomic/scripts/maf_converter.py — CIGAR creation."""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipelines", "metatranscriptomic", "scripts"))
from maf_converter import cigar_code_creation


class TestCigarCodeCreation:
    def test_all_matches(self):
        assert cigar_code_creation("ACGT", "ACGT") == "4M"

    def test_mismatch_counted_as_match(self):
        assert cigar_code_creation("ACGT", "AGGT") == "4M"

    def test_insertion_gap_in_ref(self):
        # gap in reference = insertion in read
        assert cigar_code_creation("-ACG", "TACG") == "1I3M"

    def test_deletion_gap_in_read(self):
        # gap in read = deletion from reference
        assert cigar_code_creation("ACGT", "A-GT") == "1M1D2M"

    def test_multiple_insertions(self):
        assert cigar_code_creation("--ACGT", "TTACGT") == "2I4M"

    def test_multiple_deletions(self):
        assert cigar_code_creation("ACGT", "A--T") == "1M2D1M"

    def test_single_base_match(self):
        assert cigar_code_creation("A", "A") == "1M"

    def test_empty_strings(self):
        assert cigar_code_creation("", "") == ""
