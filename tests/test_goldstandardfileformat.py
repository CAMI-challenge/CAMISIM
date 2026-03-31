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


class TestGetDictUniqueIdToGenomeFilePath:
    def test_basic_mapping(self, gsf, tmp_path):
        f = tmp_path / "mapping.tsv"
        f.write_text("genome1\t/path/to/g1.fa\ngenome2\t/path/to/g2.fa\n")
        result = gsf.get_dict_unique_id_to_genome_file_path(str(f))
        assert result == {"genome1": "/path/to/g1.fa", "genome2": "/path/to/g2.fa"}


class TestGetDictGenomeIdToTaxId:
    def test_basic_metadata(self, gsf, tmp_path):
        f = tmp_path / "metadata.tsv"
        f.write_text("genome_ID\tOTU\tNCBI_ID\tnovelty\ng1\totu1\t12345\tknown\ng2\totu2\t67890\tnovel\n")
        result = gsf.get_dict_genome_id_to_tax_id(str(f))
        assert result == {"g1": "12345", "g2": "67890"}


class TestGetDictAnonymousToOriginalId:
    def test_reverse_mapping(self, gsf, tmp_path):
        f = tmp_path / "anon_map.tsv"
        f.write_text("original_seq1\tanon_0\noriginal_seq2\tanon_1\n")
        result = gsf.get_dict_anonymous_to_original_id(str(f))
        assert result == {"anon_0": "original_seq1", "anon_1": "original_seq2"}


class TestGetDictSequenceNameToAnonymous:
    def test_forward_mapping(self, gsf, tmp_path):
        f = tmp_path / "seq_map.tsv"
        f.write_text("original_seq1\tanon_0\noriginal_seq2\tanon_1\n")
        result = gsf.get_dict_sequence_name_to_anonymous(str(f))
        assert result == {"original_seq1": "anon_0", "original_seq2": "anon_1"}


class TestGetDictSequenceNameToPositions:
    def test_basic_positions(self, gsf, tmp_path):
        f = tmp_path / "positions.tsv"
        f.write_text("seq1-0\t100\nseq1-1\t200\nseq2-0\t50\n")
        result = gsf.get_dict_sequence_name_to_positions([str(f)])
        assert "seq1" in result
        assert sorted(result["seq1"]) == [100, 200]
        assert result["seq2"] == [50]

    def test_wgsim_format(self, gsf, tmp_path):
        f = tmp_path / "positions.tsv"
        f.write_text("seq1_100_200_0:0:0_0:0:0_1\t150\n")
        result = gsf.get_dict_sequence_name_to_positions([str(f)], wgsim=True)
        assert "seq1" in result

    def test_multiple_files(self, gsf, tmp_path):
        f1 = tmp_path / "pos1.tsv"
        f1.write_text("seq1-0\t100\n")
        f2 = tmp_path / "pos2.tsv"
        f2.write_text("seq1-1\t200\n")
        result = gsf.get_dict_sequence_name_to_positions([str(f1), str(f2)])
        assert sorted(result["seq1"]) == [100, 200]
