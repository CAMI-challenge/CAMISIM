"""Tests for pipelines/metagenomic/scripts/clean_up_sequences.py"""

import pytest
from clean_up_sequences import get_new_name


class TestGetNewName:
    def test_no_collision(self):
        assert get_new_name("seq1", set()) == "seq1"

    def test_single_collision(self):
        names = {"seq1"}
        result = get_new_name("seq1", names)
        assert result == "seq1_0"

    def test_multiple_collisions(self):
        names = {"seq1", "seq1_0", "seq1_1"}
        result = get_new_name("seq1", names)
        assert result == "seq1_2"

    def test_no_collision_with_other_names(self):
        names = {"seq2", "seq3"}
        assert get_new_name("seq1", names) == "seq1"

    def test_result_not_in_existing_set(self):
        names = {"contig", "contig_0"}
        result = get_new_name("contig", names)
        assert result not in names
