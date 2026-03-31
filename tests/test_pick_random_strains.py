"""Tests for pipelines/metagenomic/scripts/pick_random_strains.py — helper utilities."""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipelines", "metagenomic", "scripts"))
from pick_random_strains import get_empty_row


class TestGetEmptyRow:
    def test_returns_dict_by_default(self):
        result = get_empty_row(["a", "b", "c"])
        assert result == {"a": "", "b": "", "c": ""}

    def test_custom_default_value(self):
        result = get_empty_row(["x", "y"], default_value="NA")
        assert result == {"x": "NA", "y": "NA"}

    def test_as_list(self):
        result = get_empty_row(["a", "b", "c"], as_list=True)
        assert result == ["", "", ""]

    def test_as_list_with_default(self):
        result = get_empty_row(["a", "b"], default_value="0", as_list=True)
        assert result == ["0", "0"]

    def test_empty_columns(self):
        assert get_empty_row([]) == {}
        assert get_empty_row([], as_list=True) == []

    def test_non_string_default_raises(self):
        with pytest.raises(AssertionError):
            get_empty_row(["a"], default_value=0)
