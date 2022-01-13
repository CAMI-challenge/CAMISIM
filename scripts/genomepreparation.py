__author__ = 'peter hofmann'
__version__ = '0.0.2'


import os
from Bio import SeqIO
from scripts.Validator.sequencevalidator import SequenceValidator
from scripts.MetaDataTable.metadatatable import MetadataTable


# ##################################
#
#          GenomePreparation
#
# ##################################


class GenomePreparation(SequenceValidator):

	_filename_seq_map = "sequence_id_map.txt"

	def __init__(self, label="GenomePreparation", logfile=None, verbose=False, debug=False):
		super(GenomePreparation, self).__init__(label=label, logfile=logfile, verbose=verbose, debug=debug)

	def write_genome_id_to_path_map(self, genome_id_to_path_map, file_path_output):
		"""
		Write mapping of genome id to genome file path to a file.

		@param file_path_output: File path
		@type file_path_output: file | FileIO | StringIO
		@param genome_id_to_path_map:
		@type genome_id_to_path_map: dict[str|unicode, str|unicode]
		"""
		with open(file_path_output, 'w') as stream_out:
			self._stream_genome_id_to_path_map(stream_out, genome_id_to_path_map)

	def _stream_genome_id_to_path_map(self, stream_out, genome_id_to_path_map):
		"""
		Write mapping of genome id to genome file path to a stream.

		@param stream_out: Stream like object
		@type stream_out: file | FileIO | StringIO
		@param genome_id_to_path_map:
		@type genome_id_to_path_map: dict[str|unicode, str|unicode]
		"""
		assert self.is_stream(stream_out)
		for genome_id, file_path in genome_id_to_path_map.items():
			stream_out.write("{}\t{}\n".format(genome_id, file_path))

	def _get_genome_id_to_path_map(self, file_path_of_file_mapping_genome_id_to_paths, list_of_drawn_genome_id):
		"""
		Get a dictionary mapping genome id to the path of their genome

		@param file_path_of_file_mapping_genome_id_to_paths: File path to file with format 'id \t path'
		@type file_path_of_file_mapping_genome_id_to_paths: str | unicode
		@param list_of_drawn_genome_id: List of genome identifiers
		@type list_of_drawn_genome_id: list[str|unicode]

		@return: genome ids mapped to their gnome file path
		@rtype: dict[str|unicode, str|unicode]
		"""
		genome_id_to_path_map = {}
		mdt = MetadataTable(logfile=self._logfile, verbose=self._verbose)
		mdt.read(file_path_of_file_mapping_genome_id_to_paths)
		if mdt.get_number_of_rows() > 0:
			genome_id_to_path_map = mdt.get_map(0, 1, unique_key=True)
		msg = "'{}' is missing one or more genome id".format(os.path.basename(file_path_of_file_mapping_genome_id_to_paths))
		assert set(genome_id_to_path_map.keys()).issuperset(list_of_drawn_genome_id), msg
		return {genome_id: genome_id_to_path_map[genome_id] for genome_id in list_of_drawn_genome_id}

	def _move_genome_file(
		self, file_path_input, file_path_output,
		stream_map, genome_id, sequence_min_length=1, set_of_sequence_names=None, file_format="fasta"):
		"""
		Move genomes into project folder, cleaning it up in the process.
		Makes sure sequence ids are unique and descriptions/comments are removed
		Returns total length of all sequences

		@attention file_format: Anything but 'fasta' is not supported, yet

		@param file_path_input: File path to raw genome
		@type file_path_input: str | unicode
		@param file_path_output: Destination path
		@type file_path_output: str | unicode
		@param sequence_min_length: Minimum length of sequences
		@type sequence_min_length: int | long
		@param set_of_sequence_names: Set of all previously used sequence names, making sure all will be unique
		@type set_of_sequence_names: set[str|unicode]
		@param file_format: 'fasta' format by default.
		@type file_format: str | unicode

		@raise Exception:
		"""
		assert self.validate_file(file_path_input)
		assert file_format == "fasta", "'{}' is not supported, yet.".format(file_format)
		#assert not self.validate_file(file_path_output, silent=True), "Overwriting files prohibited: '{}'".format(
		#	file_path_output)
		if set_of_sequence_names is None:
			set_of_sequence_names = []
		
		if (self.validate_file(file_path_output, silent=True)):
			self._logger.warning("File %s existing, skipping" % file_path_output)
			with open(file_path_input, 'r') as stream_input:
				self._add_sequences_to_map(stream_input, stream_map, genome_id, sequence_min_length, set_of_sequence_names)
			return
		with open(file_path_input, 'r') as stream_input, open(file_path_output, 'w') as stream_output:
			total_base_pairs = self._cleanup_and_filter_sequences(
				stream_input, stream_output, stream_map, genome_id, sequence_min_length, set_of_sequence_names, file_format)
		if total_base_pairs == 0:
			msg = "No valid sequences in '{}'".format(stream_input.name)
			self._logger.error(msg)
			raise Exception(msg)
	
	def _add_sequences_to_map(self, stream_input, stream_map, genome_id, sequence_min_length, set_of_sequence_names, file_format="fasta"):
		for seq_record in SeqIO.parse(stream_input, file_format):
			# remove description, else art illumina messes up sam format
			seq_record.description = ''
			if len(seq_record.seq) < sequence_min_length:
				self._logger.debug("'{}', Removing short sequence '{}', length: {}".format(
					os.path.basename(stream_input.name), seq_record.id, len(seq_record.seq)))
				continue
			if seq_record.id in set_of_sequence_names:
				new_id = self._get_new_name(seq_record.id, set_of_sequence_names)
				stream_map.write("{}\t{}\t{}\n".format(genome_id, seq_record.id, new_id))
				seq_record.id = new_id
			set_of_sequence_names.add(seq_record.id)

	def _cleanup_and_filter_sequences(
		self, stream_input, stream_output, stream_map,
		genome_id, sequence_min_length, set_of_sequence_names, file_format="fasta"):
		"""
		Rename ids that are not unique and remove sequences that are shorter than a given minimum

		@attention file_format: Anything but 'fasta' is not supported, yet

		@param stream_input: input stream of sequence file
		@type stream_input: file | FileIO | StringIO
		@param stream_output: Output stream
		@type stream_output: file | FileIO | StringIO
		@param stream_output: Output stream mapping
		@type stream_output: file | FileIO | StringIO
		@param sequence_min_length: Minimum length of sequences
		@type sequence_min_length: int | long
		@param set_of_sequence_names: Set of all previously used sequence names, making sure all will be unique
		@type set_of_sequence_names: set[str|unicode]
		@param file_format: 'fasta' format by default.
		@type file_format: str | unicode

		@return: Total length of all sequences (base pairs)
		@rtype: int | long
		"""
		total_base_pairs = 0
		for seq_record in SeqIO.parse(stream_input, file_format):
			# remove description, else art illumina messes up sam format
			seq_record.description = ''
			if len(seq_record.seq) < sequence_min_length:
				self._logger.debug("'{}', Removing short sequence '{}', length: {}".format(
					os.path.basename(stream_input.name), seq_record.id, len(seq_record.seq)))
				continue
			if seq_record.id in set_of_sequence_names:
				new_id = self._get_new_name(seq_record.id, set_of_sequence_names)
				stream_map.write("{}\t{}\t{}\n".format(genome_id, seq_record.id, new_id))
				seq_record.id = new_id
			set_of_sequence_names.add(seq_record.id)
			# file_handler.write(">{}\n".format(sequence_id))
			# file_handler.writelines("{}\n".format(seq_record.seq))
			stream_output.write(seq_record.format(file_format))
			total_base_pairs += len(seq_record.seq)
		return total_base_pairs

	def _stream_sequences_of_min_length(
		self, stream_input, stream_output, sequence_min_length, file_format="fasta"):
		"""
		Stream sequences of a minimum length

		@attention file_format: Anything but 'fasta' is not supported, yet

		@param stream_input: input stream of sequence file
		@type stream_input: file | FileIO | StringIO
		@param stream_output: Output stream
		@type stream_output: file | FileIO | StringIO
		@param stream_output: Output stream mapping
		@type stream_output: file | FileIO | StringIO
		@param sequence_min_length: Minimum length of sequences
		@type sequence_min_length: int | long
		@param file_format: 'fasta' format by default.
		@type file_format: str | unicode

		@return: Total length of all sequences (base pairs)
		@rtype: int | long
		"""
		total_base_pairs = 0
		for seq_record in SeqIO.parse(stream_input, file_format):
			# remove description, else art illumina messes up sam format
			seq_record.description = ''
			if len(seq_record.seq) < sequence_min_length:
				self._logger.debug("'{}', Removing short sequence '{}', length: {}".format(
					os.path.basename(stream_input.name), seq_record.id, len(seq_record.seq)))
				continue
			stream_output.write(seq_record.format(file_format))
			total_base_pairs += len(seq_record.seq)
		return total_base_pairs

	def move_genome_files(
		self, genome_id_to_path_map, directory_output, sequence_min_length=0, set_of_sequence_names=None):
		"""
		Move and clean up a list of genomes

		@param genome_id_to_path_map: Dictionary with file path by genome ids
		@type genome_id_to_path_map: dict[str|unicode, str|unicode]
		@param directory_output: Output directory
		@type directory_output: str | unicode
		@param sequence_min_length: Minimum length of sequences
		@type sequence_min_length: int | long
		@param set_of_sequence_names: Set of all previously used sequence names, making sure all will be unique
		@type set_of_sequence_names: set[str|unicode]
		"""
		directory_output = self.get_full_path(directory_output)
		assert isinstance(genome_id_to_path_map, dict)
		if set_of_sequence_names is None:
			set_of_sequence_names = set()
		file_path_sequence_map = os.path.join(directory_output, self._filename_seq_map)
		with open(file_path_sequence_map, 'w') as stream_map:
			for genome_id, genome_file_path in genome_id_to_path_map.items():
				file_name = os.path.basename(genome_file_path)
				new_genome_file_path = os.path.join(directory_output, file_name)
				self._move_genome_file(
					genome_file_path, new_genome_file_path, stream_map, genome_id, sequence_min_length, set_of_sequence_names)
				genome_id_to_path_map[genome_id] = new_genome_file_path

	@staticmethod
	def _get_new_name(name, set_of_sequence_names):
		"""
		Get a unique sequence name

		@param name: Current sequence name
		@type name: str | unicode
		@param set_of_sequence_names: Set of all previously used sequence names, making sure all will be unique
		@type set_of_sequence_names: set[str|unicode]
		@return: Unique sequence name
		@rtype : str | unicode
		"""
		index = 0
		new_name = name
		while new_name in set_of_sequence_names:
			new_name = name + "_{}".format(index)
			index += 1
		return new_name

	def validate_format(self, list_of_file_paths, file_format="fasta", sequence_type="dna", ambiguous=True):
		"""
		Validate file format of a list of fasta files

		@param list_of_file_paths: List of fasta file paths
		@type list_of_file_paths: list[str|unicode]
		@param file_format: 'fasta' or 'fastq'
		@type file_format: str | unicode
		@param sequence_type: 'dna' or 'rna' or 'protein'
		@type sequence_type: str | unicode
		@param ambiguous: If true ambiguous characters are valid
		@type ambiguous: bool

		@return: True if all valid
		@rtype: bool
		"""
		result = True
		for file_path in list_of_file_paths:
			if not self.validate_sequence_file(file_path, file_format, sequence_type, ambiguous):
				result = False
		return result

	def get_sequence_lengths(
		self, file_path, file_format, sequence_type, ambiguous, key=None, silent=False):
		"""
			Validate a file to be correctly formatted and have sequences of a minimum length

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
			@rtype: tuple[int|long, int|long]
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
		total_length = 0
		sequence_count = 0
		min_sequence_length = None
		with open(file_path) as file_handle:
			if not self._validate_file_start(file_handle, file_format):
				msg = "{}Invalid beginning of file '{}'.".format(prefix, os.path.basename(file_path))
				self._logger.error(msg)
				raise IOError(msg)
			try:
				for seq_record in SeqIO.parse(file_handle, file_format):
					sequence_count += 1
					if not self._validate_sequence_record(seq_record, set_of_seq_id, file_format, key=None, silent=False):
						msg = "{}{}. sequence '{}' is invalid.".format(prefix, sequence_count, seq_record.id)
						self._logger.error(msg)
						raise IOError(msg)
					if not min_sequence_length:
						min_sequence_length = len(seq_record.seq)
					total_length += len(seq_record.seq)
					if len(seq_record.seq) < min_sequence_length:
						min_sequence_length = len(seq_record.seq)
			except Exception as e:
				if not silent:
					msg = "{}Corrupt sequence in file '{}'.\nException: {}".format(prefix, os.path.basename(file_path), e.message)
					self._logger.error(msg)
					raise IOError(msg)
				return False
		if sequence_count == 0:
			return 0, 0
		return min_sequence_length, total_length
