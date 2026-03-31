"""Tests for pipelines/metagenomic/scripts/get_genomes.py — pure utility functions."""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipelines", "metagenomic", "scripts"))

pytest.importorskip("biom", reason="biom not installed")
from get_genomes import sort_by_abundance, read_genomes_list, str2bool, truncated_geometric


class TestSortByAbundance:
    def test_sorts_descending(self):
        profile = {
            "otu1": (["k__Bacteria"], [1.0, 2.0]),
            "otu2": (["k__Archaea"], [10.0, 20.0]),
            "otu3": (["k__Fungi"], [5.0, 5.0]),
        }
        result = sort_by_abundance(profile)
        assert result == ["otu2", "otu3", "otu1"]

    def test_single_otu(self):
        profile = {"otu1": ([], [3.0])}
        assert sort_by_abundance(profile) == ["otu1"]

    def test_equal_abundance(self):
        profile = {
            "a": ([], [5.0]),
            "b": ([], [5.0]),
        }
        result = sort_by_abundance(profile)
        assert len(result) == 2
        assert set(result) == {"a", "b"}


class TestReadGenomesList:
    def test_basic_parsing(self, tmp_path):
        f = tmp_path / "genomes.tsv"
        f.write_text("12345\tE. coli\tftp://example.com/path\n")
        result, total = read_genomes_list(str(f))
        assert total == 1
        assert "12345" in result
        name, paths, novelty = result["12345"]
        assert name == "E. coli"
        assert "http://example.com/path" in paths[0]  # ftp replaced with http
        assert novelty == "known_strain"

    def test_multiple_genomes_same_taxid(self, tmp_path):
        f = tmp_path / "genomes.tsv"
        f.write_text("100\tSpeciesA\tftp://a/path1\n100\tSpeciesA\tftp://a/path2\n")
        result, total = read_genomes_list(str(f))
        assert total == 2
        assert len(result["100"][1]) == 2

    def test_with_additional_file(self, tmp_path):
        main = tmp_path / "genomes.tsv"
        main.write_text("100\tSpeciesA\tftp://a/path\n")
        add = tmp_path / "additional.tsv"
        add.write_text("200\tSpeciesB\t/local/genome.fa\tnovel_strain\n")
        result, total = read_genomes_list(str(main), str(add))
        assert total == 2
        assert "200" in result
        assert result["200"][2] == "novel_strain"


class TestStr2Bool:
    @pytest.mark.parametrize("val", ["1", "true", "True", "TRUE", "yes", "y", "Y"])
    def test_truthy(self, val):
        assert str2bool(val) is True

    @pytest.mark.parametrize("val", ["0", "false", "False", "no", "n", ""])
    def test_falsy(self, val):
        assert str2bool(val) is False


class TestTruncatedGeometric:
    def test_within_bounds(self):
        for _ in range(50):
            val = truncated_geometric(0.5, 1, 5)
            assert 1 <= val <= 5

    def test_single_value_range(self):
        val = truncated_geometric(0.5, 3, 3)
        assert val == 3

    def test_invalid_params_raises(self):
        with pytest.raises(ValueError):
            truncated_geometric(0, 1, 5)

    def test_invalid_range_raises(self):
        with pytest.raises(ValueError):
            truncated_geometric(0.5, 10, 1)
