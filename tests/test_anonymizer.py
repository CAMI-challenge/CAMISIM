"""Tests for pipelines/shared/scripts/anonymizer.py — sequence anonymization."""

import io
import pytest
from anonymizer import Anonymizer


@pytest.fixture
def anon():
    return Anonymizer()


def _make_fasta(*seqs):
    """Create a FASTA-formatted StringIO from (id, seq) tuples."""
    lines = []
    for sid, seq in seqs:
        lines.append(f">{sid}\n{seq}\n")
    return io.StringIO("".join(lines))


class TestAnonymizeSequences:
    def test_basic_anonymization(self, anon):
        inp = _make_fasta(("seq1", "ACGT"), ("seq2", "TTTT"))
        out = io.StringIO()
        mapping = io.StringIO()
        anon.anonymize_sequences(mapping, input_stream=inp, output_stream=out, file_format="fasta")

        mapping.seek(0)
        lines = mapping.read().strip().split("\n")
        assert len(lines) == 2
        assert lines[0].startswith("seq1\t")
        assert lines[1].startswith("seq2\t")

    def test_prefix(self, anon):
        inp = _make_fasta(("seq1", "ACGT"))
        out = io.StringIO()
        mapping = io.StringIO()
        anon.anonymize_sequences(mapping, input_stream=inp, output_stream=out, sequence_prefix="S0R", file_format="fasta")

        mapping.seek(0)
        assert "S0R0" in mapping.read()

    def test_sequential_ids(self, anon):
        inp = _make_fasta(("a", "A"), ("b", "C"), ("c", "G"))
        out = io.StringIO()
        mapping = io.StringIO()
        anon.anonymize_sequences(mapping, input_stream=inp, output_stream=out, file_format="fasta")

        mapping.seek(0)
        ids = [line.split("\t")[1] for line in mapping.read().strip().split("\n")]
        assert ids == ["0", "1", "2"]

    def test_output_contains_sequences(self, anon):
        inp = _make_fasta(("seq1", "ACGT"))
        out = io.StringIO()
        mapping = io.StringIO()
        anon.anonymize_sequences(mapping, input_stream=inp, output_stream=out, file_format="fasta")

        out.seek(0)
        content = out.read()
        assert "ACGT" in content

    def test_invalid_format_raises(self, anon):
        with pytest.raises(AssertionError):
            anon.anonymize_sequences(io.StringIO(), file_format="bam")


class TestAnonymizeSequencePairs:
    def test_paired_ids(self, anon):
        inp = _make_fasta(("r1/1", "ACGT"), ("r1/2", "TTTT"), ("r2/1", "GGGG"), ("r2/2", "CCCC"))
        out = io.StringIO()
        mapping = io.StringIO()
        anon.anonymize_sequence_pairs(mapping, input_stream=inp, output_stream=out, file_format="fasta")

        mapping.seek(0)
        new_ids = [line.split("\t")[1] for line in mapping.read().strip().split("\n")]
        # pair 0: 0/1, 0/2; pair 1: 1/1, 1/2
        assert new_ids == ["0/1", "0/2", "1/1", "1/2"]

    def test_prefix_with_pairs(self, anon):
        inp = _make_fasta(("a", "A"), ("b", "C"))
        out = io.StringIO()
        mapping = io.StringIO()
        anon.anonymize_sequence_pairs(mapping, input_stream=inp, output_stream=out, sequence_prefix="X", file_format="fasta")

        mapping.seek(0)
        new_ids = [line.split("\t")[1] for line in mapping.read().strip().split("\n")]
        assert new_ids == ["X0/1", "X0/2"]

    def test_six_sequences_three_pairs(self, anon):
        inp = _make_fasta(
            ("r1f", "AAAA"), ("r1r", "TTTT"),
            ("r2f", "GGGG"), ("r2r", "CCCC"),
            ("r3f", "ACGT"), ("r3r", "TGCA"),
        )
        out = io.StringIO()
        mapping = io.StringIO()
        anon.anonymize_sequence_pairs(mapping, input_stream=inp, output_stream=out, file_format="fasta")

        mapping.seek(0)
        new_ids = [line.split("\t")[1] for line in mapping.read().strip().split("\n")]
        assert new_ids == ["0/1", "0/2", "1/1", "1/2", "2/1", "2/2"]


class TestAnonymizerEdgeCases:
    def test_empty_input(self, anon):
        inp = io.StringIO("")
        out = io.StringIO()
        mapping = io.StringIO()
        anon.anonymize_sequences(mapping, input_stream=inp, output_stream=out, file_format="fasta")
        assert out.getvalue() == ""
        assert mapping.getvalue() == ""

    def test_original_ids_in_mapping(self, anon):
        inp = _make_fasta(("my_original_id", "ACGT"))
        out = io.StringIO()
        mapping = io.StringIO()
        anon.anonymize_sequences(mapping, input_stream=inp, output_stream=out, file_format="fasta")
        mapping.seek(0)
        assert "my_original_id" in mapping.read()
