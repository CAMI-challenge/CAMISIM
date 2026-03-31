"""Tests for pipelines/metagenomic/scripts/build_ncbi_taxonomy.py — taxonomy utilities."""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipelines", "metagenomic", "scripts"))
import build_ncbi_taxonomy as tax


@pytest.fixture(autouse=True)
def _reset_globals():
    """Reset global state before each test."""
    tax.TAXID_TO_NAME.clear()
    tax.TAXID_TO_PARENT_TAXID.clear()
    tax.TAXID_TO_RANK.clear()
    tax.TAXID_OLD_TO_TAXID_NEW.clear()
    tax.IS_LEGACY = False
    yield
    tax.TAXID_TO_NAME.clear()
    tax.TAXID_TO_PARENT_TAXID.clear()
    tax.TAXID_TO_RANK.clear()
    tax.TAXID_OLD_TO_TAXID_NEW.clear()
    tax.IS_LEGACY = False


def _setup_simple_taxonomy():
    """Set up a minimal taxonomy: 1 -> 2 (superkingdom) -> 100 (genus) -> 200 (species)."""
    tax.TAXID_TO_PARENT_TAXID.update({
        "1": "1",
        "2": "1",
        "100": "2",
        "200": "100",
    })
    tax.TAXID_TO_RANK.update({
        "1": "no rank",
        "2": "superkingdom",
        "100": "genus",
        "200": "species",
    })
    tax.TAXID_TO_NAME.update({
        "1": "root",
        "2": "Bacteria",
        "100": "Escherichia",
        "200": "Escherichia coli",
    })


class TestGetUpdatedTaxid:
    def test_known_taxid(self):
        tax.TAXID_TO_RANK["12345"] = "species"
        assert tax.get_updated_taxid("12345") == "12345"

    def test_merged_taxid(self):
        tax.TAXID_TO_RANK["999"] = "species"
        tax.TAXID_OLD_TO_TAXID_NEW["111"] = "999"
        assert tax.get_updated_taxid("111") == "999"

    def test_special_taxid_other_entries(self):
        assert tax.get_updated_taxid(tax.OTHER_ENTRIES_TAXID) == tax.OTHER_ENTRIES_TAXID

    def test_special_taxid_plasmid(self):
        assert tax.get_updated_taxid(tax.PLASMID_SPECIES_TAXID) == tax.PLASMID_SPECIES_TAXID

    def test_unknown_taxid_raises(self):
        with pytest.raises(ValueError):
            tax.get_updated_taxid("9999999")


class TestGetScientificName:
    def test_known_name(self):
        tax.TAXID_TO_RANK["200"] = "species"
        tax.TAXID_TO_NAME["200"] = "Escherichia coli"
        assert tax.get_scientific_name("200") == "Escherichia coli"

    def test_merged_taxid_lookup(self):
        tax.TAXID_TO_RANK["200"] = "species"
        tax.TAXID_TO_NAME["200"] = "Escherichia coli"
        tax.TAXID_OLD_TO_TAXID_NEW["111"] = "200"
        assert tax.get_scientific_name("111") == "Escherichia coli"

    def test_unknown_raises(self):
        with pytest.raises(ValueError):
            tax.get_scientific_name("9999999")


class TestIsDescendantOf:
    def test_direct_descendant(self):
        _setup_simple_taxonomy()
        assert tax.is_descendant_of("200", "100") is True

    def test_transitive_descendant(self):
        _setup_simple_taxonomy()
        assert tax.is_descendant_of("200", "2") is True

    def test_not_descendant(self):
        _setup_simple_taxonomy()
        assert tax.is_descendant_of("2", "200") is False

    def test_self_is_descendant(self):
        _setup_simple_taxonomy()
        assert tax.is_descendant_of("200", "200") is True

    def test_unknown_taxid(self):
        assert tax.is_descendant_of("999", "1") is False


class TestGetOtherEntriesLineage:
    def test_basic_lineage(self):
        lineage = tax.get_other_entries_lineage("45202")
        assert lineage[0] == tax.OTHER_ENTRIES_TAXID
        assert lineage[8] == "45202"
        assert lineage[9] == ""
        assert len(lineage) == len(tax.OTHER_RANKS)

    def test_with_strain(self):
        lineage = tax.get_other_entries_lineage("45202", taxid_if_strain="99999")
        assert lineage[9] == "99999"

    def test_empty_middle_ranks(self):
        lineage = tax.get_other_entries_lineage("45202")
        for i in range(1, 8):
            assert lineage[i] == ""


class TestGetLineageOfLegalRanks:
    def test_legacy_ranks(self):
        tax.IS_LEGACY = True
        _setup_simple_taxonomy()
        lineage = tax.get_lineage_of_legal_ranks("200")
        assert len(lineage) == len(tax.LEGACY_RANKS)
        # species should be filled
        species_idx = tax.LEGACY_RANKS.index("species")
        assert lineage[species_idx] == "200"
        # genus should be filled
        genus_idx = tax.LEGACY_RANKS.index("genus")
        assert lineage[genus_idx] == "100"

    def test_missing_ranks_get_default(self):
        tax.IS_LEGACY = True
        _setup_simple_taxonomy()
        lineage = tax.get_lineage_of_legal_ranks("200", default_value="")
        # phylum etc. not in our mini taxonomy → default
        phylum_idx = tax.LEGACY_RANKS.index("phylum")
        assert lineage[phylum_idx] == ""


class TestParseStream:
    def test_basic_tsv(self):
        lines = ["a\tb\tc\n", "d\te\tf\n"]
        rows = list(tax.parse_stream(lines))
        assert rows == [["a", "b", "c"], ["d", "e", "f"]]

    def test_skips_comments(self):
        lines = ["# comment\n", "a\tb\n"]
        rows = list(tax.parse_stream(lines))
        assert len(rows) == 1

    def test_skips_empty_lines(self):
        lines = ["\n", "a\tb\n", "\n"]
        rows = list(tax.parse_stream(lines))
        assert len(rows) == 1

    def test_strips_newlines(self):
        lines = ["a\tb\n"]
        rows = list(tax.parse_stream(lines))
        assert rows[0] == ["a", "b"]
