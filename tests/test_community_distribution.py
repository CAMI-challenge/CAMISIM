"""Tests for pipelines/shared/scripts/get_community_distribution.py"""

import pytest
import get_community_distribution as cd


class TestLtZero:
    def test_positive_value_unchanged(self):
        assert cd.lt_zero(5.0) == 5.0

    def test_zero_returns_floor(self):
        assert cd.lt_zero(0) == 0.001

    def test_negative_returns_floor(self):
        assert cd.lt_zero(-3.5) == 0.001

    def test_small_positive_unchanged(self):
        assert cd.lt_zero(0.0001) == 0.0001


class TestIsBooleanState:
    @pytest.mark.parametrize("word", ["yes", "true", "on", "y", "t"])
    def test_truthy_words(self, word):
        assert cd.is_boolean_state(word) is True

    @pytest.mark.parametrize("word", ["no", "false", "off", "n", "f"])
    def test_falsy_words(self, word):
        assert cd.is_boolean_state(word) is True

    @pytest.mark.parametrize("word", ["maybe", "1", "0", "YES", "True", ""])
    def test_non_boolean_words(self, word):
        assert cd.is_boolean_state(word) is False


class TestGetBooleanState:
    @pytest.mark.parametrize("word,expected", [
        ("yes", True), ("true", True), ("on", True), ("y", True), ("t", True),
        ("no", False), ("false", False), ("off", False), ("n", False), ("f", False),
    ])
    def test_valid_words(self, word, expected):
        assert cd.get_boolean_state(word) is expected

    def test_invalid_word_raises(self):
        with pytest.raises(AssertionError):
            cd.get_boolean_state("maybe")


class TestStrToBool:
    def test_true(self):
        assert cd.str_to_bool("True") is True

    def test_false(self):
        assert cd.str_to_bool("False") is False

    def test_invalid_raises(self):
        with pytest.raises(ValueError):
            cd.str_to_bool("true")

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            cd.str_to_bool("")


class TestValidateNumber:
    def test_valid_number_no_bounds(self):
        assert cd.validate_number(5) is True

    def test_below_minimum(self):
        assert cd.validate_number(1, minimum=5) is False

    def test_above_maximum(self):
        assert cd.validate_number(10, maximum=5) is False

    def test_within_range(self):
        assert cd.validate_number(5, minimum=1, maximum=10) is True

    def test_zero_excluded(self):
        assert cd.validate_number(0, zero=False) is False

    def test_zero_allowed_by_default(self):
        assert cd.validate_number(0) is True

    def test_float_valid(self):
        assert cd.validate_number(3.14, minimum=1.0, maximum=4.0) is True

    def test_non_number_raises(self):
        with pytest.raises(AssertionError):
            cd.validate_number("five")


class TestGetInitialList:
    def test_dimensions(self):
        result = cd.get_initial_list(3, 2)
        assert len(result) == 3
        assert all(len(row) == 2 for row in result)

    def test_all_zeros(self):
        result = cd.get_initial_list(2, 3)
        for row in result:
            assert all(v == 0.0 for v in row)

    def test_independent_rows(self):
        result = cd.get_initial_list(2, 2)
        result[0][0] = 1.0
        assert result[1][0] == 0.0

    def test_non_int_raises(self):
        with pytest.raises(AssertionError):
            cd.get_initial_list(2.0, 3)


class TestHasUniqueColumns:
    def test_unique(self):
        assert cd.has_unique_columns(["a", "b", "c"]) is True

    def test_duplicates(self):
        assert cd.has_unique_columns(["a", "b", "a"]) is False

    def test_empty(self):
        assert cd.has_unique_columns([]) is True


class TestGetDistributionFilePaths:
    def test_correct_count(self):
        result = cd.get_distribution_file_paths("/out", 3)
        assert len(result) == 3

    def test_correct_format(self):
        result = cd.get_distribution_file_paths("/out", 2)
        assert result[0].endswith("distribution_0.txt")
        assert result[1].endswith("distribution_1.txt")

    def test_directory_prefix(self):
        result = cd.get_distribution_file_paths("/my/dir", 1)
        assert result[0].startswith("/my/dir")


class TestHasColumn:
    def test_existing_column(self):
        table = {"name": ["a", "b"], "id": [1, 2]}
        assert cd.has_column(table, "name") is True

    def test_missing_column(self):
        table = {"name": ["a", "b"]}
        assert cd.has_column(table, "id") is False

    def test_int_key(self):
        table = {0: ["a"], 1: ["b"]}
        assert cd.has_column(table, 0) is True


class TestGetMap:
    def test_basic_mapping(self):
        table = {"key": ["a", "b", "c"], "val": [1, 2, 3]}
        result = cd.get_map(table, "key", "val")
        assert result == {"a": 1, "b": 2, "c": 3}

    def test_duplicate_key_raises(self):
        table = {"key": ["a", "a"], "val": [1, 2]}
        with pytest.raises(KeyError):
            cd.get_map(table, "key", "val", unique_key=True)

    def test_missing_column_raises(self):
        table = {"key": ["a"]}
        with pytest.raises(AssertionError):
            cd.get_map(table, "key", "val")

    def test_empty_table(self):
        table = {"key": []}
        result = cd.get_map(table, "key", "key")
        assert result == {}


class TestRandomDistributionToRelativeAbundance:
    def test_sums_to_one(self):
        pop = [[3.0, 6.0], [7.0, 4.0]]
        cd.random_distribution_to_relative_abundance(pop)
        for sample_idx in range(2):
            total = sum(row[sample_idx] for row in pop)
            assert abs(total - 1.0) < 1e-9

    def test_proportions_correct(self):
        pop = [[2.0], [8.0]]
        cd.random_distribution_to_relative_abundance(pop)
        assert abs(pop[0][0] - 0.2) < 1e-9
        assert abs(pop[1][0] - 0.8) < 1e-9
