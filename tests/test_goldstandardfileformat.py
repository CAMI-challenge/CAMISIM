"""Tests for pipelines/shared/scripts/goldstandardfileformat.py — TSV parsing."""

import io
import pytest
from goldstandardfileformat import GoldStandardFileFormat


@pytest.fixture
def gsf():
    return GoldStandardFileFormat()


class TestParseColumnNames:
    def test_basic_header(self, gsf):
        stream = io.StringIO("name\tid\tvalue\nrow1\n")
        cols = gsf._parse_column_names(stream, "\t")
        assert cols == ["name", "id", "value"]

    def test_duplicate_columns_raises(self, gsf):
        stream = io.StringIO("name\tname\n")
        with pytest.raises(AssertionError):
            gsf._parse_column_names(stream, "\t")


class TestProcessLines:
    def test_without_column_names(self, gsf):
        stream = io.StringIO("a\tb\nc\td\n")
        result = gsf.process_lines(stream, {}, ["#"], "\t", column_names=False)
        assert result[0] == ["a", "c"]
        assert result[1] == ["b", "d"]

    def test_with_column_names(self, gsf):
        stream = io.StringIO("name\tid\nalice\t1\nbob\t2\n")
        result = gsf.process_lines(stream, {}, ["#"], "\t", column_names=True)
        assert result["name"] == ["alice", "bob"]
        assert result["id"] == ["1", "2"]

    def test_skips_comments(self, gsf):
        stream = io.StringIO("# comment\na\tb\n")
        result = gsf.process_lines(stream, {}, ["#"], "\t", column_names=False)
        assert result[0] == ["a"]

    def test_skips_empty_lines(self, gsf):
        stream = io.StringIO("a\tb\n\nc\td\n")
        result = gsf.process_lines(stream, {}, ["#"], "\t", column_names=False)
        assert result[0] == ["a", "c"]

    def test_bad_column_count_raises(self, gsf):
        stream = io.StringIO("name\tid\nalice\t1\nbob\n")
        with pytest.raises(ValueError, match="Bad number of values"):
            gsf.process_lines(stream, {}, ["#"], "\t", column_names=True)


class TestRead:
    def test_read_from_stream(self, gsf):
        stream = io.TextIOWrapper(io.BytesIO(b"x\ty\n1\t2\n3\t4\n"))
        result = gsf.read(stream, separator="\t", column_names=True)
        assert result["x"] == ["1", "3"]
        assert result["y"] == ["2", "4"]

    def test_read_without_headers(self, gsf, tmp_path):
        f = tmp_path / "data.tsv"
        f.write_text("a\tb\nc\td\n")
        result = gsf.read(str(f), separator="\t")
        assert result[0] == ["a", "c"]


class TestHasColumn:
    def test_existing(self, gsf):
        assert gsf.has_column({"a": [1]}, "a") is True

    def test_missing(self, gsf):
        assert gsf.has_column({"a": [1]}, "b") is False

    def test_int_key(self, gsf):
        assert gsf.has_column({0: ["x"]}, 0) is True


class TestGetMap:
    def test_basic(self, gsf):
        table = {"k": ["a", "b"], "v": [1, 2]}
        assert gsf.get_map(table, "k", "v") == {"a": 1, "b": 2}

    def test_duplicate_key_raises(self, gsf):
        table = {"k": ["a", "a"], "v": [1, 2]}
        with pytest.raises(KeyError):
            gsf.get_map(table, "k", "v", unique_key=True)

    def test_missing_column_raises(self, gsf):
        with pytest.raises(AssertionError):
            gsf.get_map({"k": []}, "k", "missing")


class TestWriteGsReadMapping:
    def test_basic_mapping(self, gsf):
        out = io.StringIO()
        anon_to_read = {"A0": "seq1.1-0", "A1": "seq1.1-1"}
        seq_to_genome = {"seq1.1": "genome1"}
        genome_to_tax = {"genome1": "12345"}
        gsf.write_gs_read_mapping(
            out, anon_to_read, seq_to_genome, genome_to_tax,
            nanosim_real_fastq=False, wgsim=False,
        )
        content = out.getvalue()
        assert "#anonymous_read_id" in content
        assert "A0\tgenome1\t12345\tseq1.1-0" in content

    def test_wgsim_mode(self, gsf):
        out = io.StringIO()
        anon_to_read = {"A0": "seq1_100_200_0:0:0_0:0:0_1/1"}
        seq_to_genome = {"seq1": "g1"}
        genome_to_tax = {"g1": "999"}
        gsf.write_gs_read_mapping(
            out, anon_to_read, seq_to_genome, genome_to_tax,
            nanosim_real_fastq=False, wgsim=True,
        )
        content = out.getvalue()
        assert "g1\t999" in content

    def test_missing_dash_raises(self, gsf):
        out = io.StringIO()
        anon_to_read = {"A0": "readwithoutdash"}
        seq_to_genome = {}
        genome_to_tax = {}
        with pytest.raises(ValueError, match="missing '-'"):
            gsf.write_gs_read_mapping(
                out, anon_to_read, seq_to_genome, genome_to_tax,
                nanosim_real_fastq=False, wgsim=False,
            )

    def test_missing_underscore_wgsim_raises(self, gsf):
        out = io.StringIO()
        anon_to_read = {"A0": "readwithoutunderscore"}
        seq_to_genome = {}
        genome_to_tax = {}
        with pytest.raises(ValueError, match="missing '_'"):
            gsf.write_gs_read_mapping(
                out, anon_to_read, seq_to_genome, genome_to_tax,
                nanosim_real_fastq=False, wgsim=True,
            )


class TestWriteGsaContigMapping:
    def test_basic_contig_mapping(self, gsf):
        out = io.StringIO()
        seq_to_anon = {"seq1_from_100_to_200_stuff": "contig_0"}
        seq_positions = {"seq1": [50, 100, 150, 200, 250]}
        seq_to_genome = {"seq1": "g1"}
        genome_to_tax = {"g1": "555"}
        gsf.write_gsa_contig_mapping(
            out, seq_to_anon, seq_positions, seq_to_genome, genome_to_tax,
        )
        content = out.getvalue()
        assert "#anonymous_contig_id" in content
        assert "contig_0\tg1\t555\tseq1" in content

    def test_bad_sequence_id_raises(self, gsf):
        out = io.StringIO()
        seq_to_anon = {"unknown_from_0_to_100_x": "c0"}
        seq_positions = {"unknown": [50]}
        seq_to_genome = {}  # missing
        genome_to_tax = {}
        with pytest.raises(KeyError):
            gsf.write_gsa_contig_mapping(
                out, seq_to_anon, seq_positions, seq_to_genome, genome_to_tax,
            )
