__author__ = 'hofmann'
__version__ = '0.0.6'

import io
import os
import string
from Bio import SeqIO
from Bio.Seq import Seq
from .validator import Validator

# Todo: allow for multiple sequence_type


class SequenceValidator(Validator):
	"""
		Colection of methods for the validation of sequences and sequence files
	"""
	_label = "SequenceValidator"

	_formats = ["fasta", "fastq"]

	_qformats = {
		"sanger": [0, 40, 33],
		"solexa": [-5, 40, 64],
		"illumina": [0, 41, 33]  # 1.8+
	}

	_sequence_indicators = {
		"fasta": ">",
		"fastq": "@"
		}

	_alphabets = {
		"rna": ["GAUC", "GAUCRYWSMKHBVDN"],
		"dna": ["GATC", "GATCRYWSMKHBVDN", "GATCRYWSMKHBVDN"],
		"protein": ["ACDEFGHIKLMNPQRSTVWYBXZJUO", "ACDEFGHIKLMNPQRSTVWYBXZJUO"]
		}

	_legal_text_characters = string.printable

	def validate_folder_with_sequence_files(
		self, directory, file_format, sequence_type, ambiguous, file_extension, key=None, silent=False):
		"""
			Validate a file to be correctly formatted

			@attention: Currently only phred quality for fastq files

			@param directory: Path to directory with files containing sequences
			@type directory: str | unicode
			@param file_format: Format of the file at the file_path provided. Valid: 'fasta', 'fastq'
			@type file_format: str | unicode
			@param sequence_type: Are the sequences DNA or RNA? Valid: 'rna', 'dna', 'protein'
			@type sequence_type: str | unicode
			@param ambiguous: True or False, DNA example for strict 'GATC',  ambiguous example 'GATCRYWSMKHBVDN'
			@type ambiguous: bool
			@param file_extension: file extension to be filtered for. Example: '.fasta' '.fq'
			@type file_extension: str | None
			@param key: If True, no error message will be made
			@type key: str | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if the file is correctly formatted
			@rtype: bool
		"""
		list_of_file_paths = self.get_files_in_directory(directory, file_extension)
		result = True
		for file_path in list_of_file_paths:
			if not self.validate_sequence_file(file_path, file_format, sequence_type, ambiguous, key, silent):
				result = False
		return result

	def _validate_sequence_record(self, seq_record, set_of_seq_id, file_format, key=None, silent=False):
			result = True
			if not self.validate_sequence(seq_record.seq, key=key, silent=silent):
				result = False
			if not self.validate_sequence_id(seq_record.id, used_ids=set_of_seq_id, key=key, silent=silent):
				result = False
			set_of_seq_id.add(seq_record.id)
			if not self.validate_sequence_description(seq_record.description, key=key, silent=silent):
				result = False
			if file_format == "fastq":
				if not self.validate_sequence_quality(
					seq_record.letter_annotations["phred_quality"], key=key, silent=silent):
					result = False
			return result

	def validate_sequence_file(self, file_path, file_format, sequence_type, ambiguous, key=None, silent=False):
		"""
			Validate a file to be correctly formatted

			@attention: Currently only phred quality for fastq files

			@param file_path: Path to file containing sequences
			@type file_path: str | unicode
			@param file_format: Format of the file at the file_path provided. Valid: 'fasta', 'fastq'
			@type file_format: str | unicode
			@param sequence_type: Are the sequences DNA or RNA? Valid: 'rna', 'dna', 'protein'
			@type sequence_type: str | unicode
			@param ambiguous: True or False, DNA example for strict 'GATC',  ambiguous example 'GATCRYWSMKHBVDN'
			@type ambiguous: bool
			@param key: If True, no error message will be made
			@type key: str | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if the file is correctly formatted
			@rtype: bool
		"""
		assert self.validate_file(file_path)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in self._formats
		assert isinstance(sequence_type, str)
		sequence_type = sequence_type.lower()
		assert sequence_type in self._alphabets

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		if ambiguous:
			alphabet = self._alphabets[sequence_type][1]
		else:
			alphabet = self._alphabets[sequence_type][0]

		set_of_seq_id = set()

		with open(file_path) as file_handle:
			if not self._validate_file_start(file_handle, file_format):
				if not silent:
					self._logger.error("{}Invalid beginning of file '{}'.".format(prefix, os.path.basename(file_path)))
				return False
			sequence_count = 0
			try:
				for seq_record in SeqIO.parse(file_handle, file_format):
					sequence_count += 1
					if not self._validate_sequence_record(seq_record, set_of_seq_id, file_format, key=None, silent=False):
						if not silent:
							self._logger.error("{}{}. sequence '{}' is invalid.".format(prefix, sequence_count, seq_record.id))
						return False
			except Exception as e:
				if not silent:
					self._logger.error("{}Corrupt sequence in file '{}'.\nException: {}".format(
						prefix, os.path.basename(file_path), e.message))
				return False
		return True

	def _validate_file_start(self, stream_handle, file_format):
		"""
			Validate that a stream with sequences starts with the correct character

			@param stream_handle: Any kind of stream type
			@type stream_handle: file | io.FileIO | StringIO.StringIO
			@param file_format: Format of the file at the file_path provided. Valid: 'fasta', 'fastq'
			@type file_format: str | unicode

			@return: True if the first character is correct
			@rtype: bool
		"""
		assert self.is_stream(stream_handle)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in self._formats

		sequence_indicator = self._sequence_indicators[file_format]

		head = stream_handle.read(1)
		stream_handle.seek(0)
		if not head:
			return False
		if not head.startswith(sequence_indicator):
			return False
		return True

	def validate_sequence_id(self, identifier, used_ids=None, key=None, silent=False):
		"""
			Validate that the sequence identifier has only valid characters

			@attention:

			@param identifier: sequence
			@type identifier: str | unicode
			@param used_ids: Set of used up ids, that should not be repeated
			@type used_ids: set
			@param key: If True, no error message will be made
			@type key: str | unicode | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if valid
			@rtype: bool
		"""
		set_of_seq_id = used_ids
		if set_of_seq_id is None:
			set_of_seq_id = set()
		assert isinstance(set_of_seq_id, set)
		assert isinstance(identifier, str)
		assert isinstance(silent, bool)

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		if not len(identifier) > 0:
			if not silent:
				self._logger.error("{}Missing sequence id".format(prefix))
			return False

		if not self.validate_characters(identifier, legal_alphabet=self._legal_text_characters, key=key, silent=silent):
			return False

		if identifier in set_of_seq_id:
			if not silent:
				self._logger.error("{}Repeated sequence id '{}'".format(
					prefix, identifier))
			return False
		return True

	def validate_sequence_description(self, description, key=None, silent=False):
		"""
			Validate that the sequence description has only valid characters

			@attention:

			@param description: sequence
			@type description: str | unicode
			@param key: If True, no error message will be made
			@type key: str | unicode | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if valid
			@rtype: bool
		"""
		assert isinstance(description, str)
		assert key is None or isinstance(key, str)
		assert isinstance(silent, bool)

		if not self.validate_characters(description, legal_alphabet=self._legal_text_characters, key=key, silent=silent):
			return False
		return True

	def validate_sequence_quality(self, quality, qformat="Illumina", key=None, silent=False):
		"""
			Validate that the sequence description has only valid characters

			@attention:

			@param quality: quality of each letter
			@type quality: list[int]
			@param qformat: 'Illumina', 'Sanger', 'Solexa'
			@type qformat: str | unicode
			@param key: If True, no error message will be made
			@type key: str | unicode | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if valid, else False
			@rtype: bool
		"""
		assert isinstance(quality, list)
		assert isinstance(silent, bool)
		assert isinstance(qformat, str)
		qformat = qformat.lower()
		assert qformat in self._qformats, "{} not in {}".format(qformat, self._qformats.keys())

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		# offset = self._qformats[qformat][2]
		minimum = self._qformats[qformat][0]  # +offset
		maximum = self._qformats[qformat][1]  # +offset

		invalid_indexes = [
			"{}: '{}'".format(index, value) for index, value in enumerate(quality) if not minimum <= value <= maximum]
		if len(invalid_indexes) > 0:
			if not silent:
				self._logger.error("{}Invalid quality at: {}.".format(prefix, ", ".join(invalid_indexes)))
			return False
		return True

	def validate_sequence(self, sequence, key=None, silent=False):
		"""
			Validate that the sequence has only valid characters

			@attention:

			@param sequence: sequence
			@type sequence: Seq
			@param key: If True, no error message will be made
			@type key: str | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if valid
			@rtype: bool
		"""
		assert isinstance(sequence, Seq)
		assert isinstance(silent, bool)

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		if not len(sequence) > 0:
			if not silent:
				self._logger.error("{}Empty sequence".format(prefix))
			return False

		if not self.validate_characters(
			sequence.upper(), key=key, silent=silent):
			return False
		return True
