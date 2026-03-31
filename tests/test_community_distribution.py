"""Tests for pipelines/shared/scripts/get_community_distribution.py"""

import io
import random
import pytest
import get_community_distribution as cd


# ---------------------------------------------------------------------------
# lt_zero
# ---------------------------------------------------------------------------
class TestLtZero:
    def test_positive_value_unchanged(self):
        assert cd.lt_zero(5.0) == 5.0

    def test_zero_returns_floor(self):
        assert cd.lt_zero(0) == 0.001

    def test_negative_returns_floor(self):
        assert cd.lt_zero(-3.5) == 0.001

    def test_small_positive_unchanged(self):
        assert cd.lt_zero(0.0001) == 0.0001


# ---------------------------------------------------------------------------
# Boolean helpers
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# validate_number
# ---------------------------------------------------------------------------
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

    def test_at_exact_minimum(self):
        assert cd.validate_number(5, minimum=5) is True

    def test_at_exact_maximum(self):
        assert cd.validate_number(5, maximum=5) is True

    def test_negative_number(self):
        assert cd.validate_number(-3, minimum=-10, maximum=0) is True


# ---------------------------------------------------------------------------
# get_initial_list
# ---------------------------------------------------------------------------
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

    def test_single_element(self):
        result = cd.get_initial_list(1, 1)
        assert result == [[0.0]]


# ---------------------------------------------------------------------------
# has_unique_columns
# ---------------------------------------------------------------------------
class TestHasUniqueColumns:
    def test_unique(self):
        assert cd.has_unique_columns(["a", "b", "c"]) is True

    def test_duplicates(self):
        assert cd.has_unique_columns(["a", "b", "a"]) is False

    def test_empty(self):
        assert cd.has_unique_columns([]) is True

    def test_single(self):
        assert cd.has_unique_columns(["a"]) is True


# ---------------------------------------------------------------------------
# get_distribution_file_paths
# ---------------------------------------------------------------------------
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

    def test_zero_samples(self):
        result = cd.get_distribution_file_paths("/out", 0)
        assert result == []


# ---------------------------------------------------------------------------
# get_compression_type
# ---------------------------------------------------------------------------
class TestGetCompressionType:
    def test_gz(self):
        assert cd.get_compression_type("data.fastq.gz") == "gz"

    def test_bz2(self):
        assert cd.get_compression_type("data.tar.bz2") == "bz2"

    def test_no_compression(self):
        assert cd.get_compression_type("data.txt") is None

    def test_unknown_extension(self):
        assert cd.get_compression_type("data.xyz") is None


# ---------------------------------------------------------------------------
# has_column / get_map
# ---------------------------------------------------------------------------
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

    def test_same_column_for_key_and_value(self):
        # get_map requires len(metadata_table) >= 2 columns to produce results
        table = {"col": ["x", "y", "z"], "other": ["a", "b", "c"]}
        result = cd.get_map(table, "col", "col")
        assert result == {"x": "x", "y": "y", "z": "z"}


# ---------------------------------------------------------------------------
# Distribution functions (seeded for reproducibility)
# ---------------------------------------------------------------------------
class TestAddInitialLogDistribution:
    def test_fills_first_sample(self):
        pop = cd.get_initial_list(5, 3)
        random.seed(42)
        cd.add_initial_log_distribution(pop, 1.0, 1.0)
        for row in pop:
            assert row[0] > 0  # lognormal is always positive
            assert row[1] == 0.0  # other samples untouched
            assert row[2] == 0.0

    def test_reproducible_with_seed(self):
        pop1 = cd.get_initial_list(3, 1)
        random.seed(99)
        cd.add_initial_log_distribution(pop1, 2.0, 0.5)

        pop2 = cd.get_initial_list(3, 1)
        random.seed(99)
        cd.add_initial_log_distribution(pop2, 2.0, 0.5)
        assert pop1 == pop2


class TestAddReplicates:
    def test_fills_remaining_samples(self):
        pop = cd.get_initial_list(3, 4)
        random.seed(42)
        cd.add_initial_log_distribution(pop, 1.0, 1.0)
        cd.add_replicates(pop, 0, 0.1)
        for row in pop:
            assert all(v > 0 for v in row)

    def test_first_sample_unchanged(self):
        pop = cd.get_initial_list(2, 3)
        random.seed(42)
        cd.add_initial_log_distribution(pop, 1.0, 1.0)
        first_vals = [row[0] for row in pop]
        cd.add_replicates(pop, 0, 0.1)
        for i, row in enumerate(pop):
            assert row[0] == first_vals[i]


class TestAddTimeseriesGauss:
    def test_extinction_propagates(self):
        pop = [[0.0, 0.0, 0.0]]
        cd.add_timeseries_gauss(pop, 0, 1.0)
        assert pop[0][1] == 0.0
        assert pop[0][2] == 0.0

    def test_positive_values_produce_positive(self):
        pop = [[5.0, 0.0, 0.0]]
        random.seed(42)
        cd.add_timeseries_gauss(pop, 0, 0.001)
        assert all(v > 0 for v in pop[0])


class TestAddTimeseriesLognorm:
    def test_fills_sequentially(self):
        pop = cd.get_initial_list(2, 3)
        random.seed(42)
        cd.add_initial_log_distribution(pop, 1.0, 1.0)
        cd.add_timeseries_lognorm(pop, 1.0, 1.0)
        for row in pop:
            assert all(v > 0 for v in row)


class TestAddDifferential:
    def test_fills_independently(self):
        pop = cd.get_initial_list(2, 3)
        random.seed(42)
        cd.add_initial_log_distribution(pop, 1.0, 1.0)
        cd.add_differential(pop, 1.0, 1.0)
        for row in pop:
            assert all(v > 0 for v in row[1:])


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

    def test_single_genome(self):
        pop = [[5.0, 3.0]]
        cd.random_distribution_to_relative_abundance(pop)
        assert pop[0][0] == 1.0
        assert pop[0][1] == 1.0


# ---------------------------------------------------------------------------
# get_lists_of_distributions (integration of the above)
# ---------------------------------------------------------------------------
class TestGetListsOfDistributions:
    @pytest.mark.parametrize("modus", [
        "replicates", "timeseries_normal", "timeseries_lognormal", "differential",
    ])
    def test_all_modes_produce_valid_distributions(self, modus):
        random.seed(42)
        result = cd.get_lists_of_distributions(
            size_of_population=5,
            number_of_samples=3,
            modus=modus,
            log_mu=1.0,
            log_sigma=1.0,
            gauss_mu=0.0,
            gauss_sigma=0.5,
        )
        assert len(result) == 5
        assert all(len(row) == 3 for row in result)
        # each sample column should sum to ~1.0
        for s in range(3):
            total = sum(row[s] for row in result)
            assert abs(total - 1.0) < 1e-6

    def test_deterministic_with_seed(self):
        random.seed(123)
        r1 = cd.get_lists_of_distributions(3, 2, "replicates", 1.0, 1.0, 0.0, 0.5)
        random.seed(123)
        r2 = cd.get_lists_of_distributions(3, 2, "replicates", 1.0, 1.0, 0.0, 0.5)
        assert r1 == r2


# ---------------------------------------------------------------------------
# write_distribution_file
# ---------------------------------------------------------------------------
class TestWriteDistributionFile:
    def test_basic_output(self):
        out = io.StringIO()
        abundances = {"genome_A": [0.3, 0.7], "genome_B": [0.6, 0.4]}
        cd.write_distribution_file(out, abundances, 1)
        content = out.getvalue()
        assert "genome_A\t0.3\n" in content
        assert "genome_B\t0.6\n" in content

    def test_second_sample(self):
        out = io.StringIO()
        abundances = {"g1": [0.1, 0.9]}
        cd.write_distribution_file(out, abundances, 2)
        assert "g1\t0.9\n" in out.getvalue()

    def test_all_genomes_written(self):
        out = io.StringIO()
        abundances = {f"g{i}": [float(i)] for i in range(10)}
        cd.write_distribution_file(out, abundances, 1)
        lines = [l for l in out.getvalue().strip().split("\n") if l]
        assert len(lines) == 10


# ---------------------------------------------------------------------------
# read (file-based TSV reader)
# ---------------------------------------------------------------------------
class TestRead:
    def test_basic_tsv(self, tmp_path):
        f = tmp_path / "data.tsv"
        f.write_text("a\t1\nb\t2\nc\t3\n")
        result = cd.read(str(f))
        assert result[0] == ["a", "b", "c"]
        assert result[1] == ["1", "2", "3"]

    def test_skips_comments(self, tmp_path):
        f = tmp_path / "data.tsv"
        f.write_text("# header\na\t1\nb\t2\n")
        result = cd.read(str(f))
        assert result[0] == ["a", "b"]

    def test_skips_empty_lines(self, tmp_path):
        f = tmp_path / "data.tsv"
        f.write_text("a\t1\n\nb\t2\n")
        result = cd.read(str(f))
        assert result[0] == ["a", "b"]

    def test_custom_separator(self, tmp_path):
        f = tmp_path / "data.csv"
        f.write_text("a,1\nb,2\n")
        result = cd.read(str(f), separator=",")
        assert result[0] == ["a", "b"]
        assert result[1] == ["1", "2"]

    def test_custom_comment_char(self, tmp_path):
        f = tmp_path / "data.tsv"
        f.write_text("% comment\na\t1\n")
        result = cd.read(str(f), comment_line="%")
        assert result[0] == ["a"]


# ---------------------------------------------------------------------------
# openy (file opener with compression support)
# ---------------------------------------------------------------------------
class TestOpeny:
    def test_open_plain_text(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("hello")
        fh = cd.openy(str(f))
        assert fh.read() == "hello"
        fh.close()

    def test_invalid_mode_raises(self):
        with pytest.raises(AssertionError):
            cd.openy("test.txt", mode="x")

    def test_write_and_read_gz(self, tmp_path):
        f = str(tmp_path / "test.gz")
        fh = cd.openy(f, mode="w", compression_type="gz")
        fh.write(b"compressed data")
        fh.close()
        fh = cd.openy(f, mode="r")
        assert b"compressed data" in fh.read()
        fh.close()


# ---------------------------------------------------------------------------
# get_genome_id_to_path_map
# ---------------------------------------------------------------------------
class TestGetGenomeIdToPathMap:
    def test_basic_mapping(self, tmp_path):
        f = tmp_path / "genome_locations.tsv"
        f.write_text("genome1\t/path/to/g1.fa\ngenome2\t/path/to/g2.fa\n")
        result = cd.get_genome_id_to_path_map(str(f))
        assert result == {"genome1": "/path/to/g1.fa", "genome2": "/path/to/g2.fa"}

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.tsv"
        f.write_text("")
        result = cd.get_genome_id_to_path_map(str(f))
        assert result == {}

    def test_single_genome(self, tmp_path):
        f = tmp_path / "single.tsv"
        f.write_text("gA\t/genomes/a.fasta\n")
        result = cd.get_genome_id_to_path_map(str(f))
        assert result == {"gA": "/genomes/a.fasta"}
