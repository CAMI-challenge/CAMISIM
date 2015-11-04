__author__ = 'Peter Hofmann'

import os
import time
import shutil
import tempfile
import scripts.parallel as parallel
from scripts.MGExtract.sequencemerger import SequenceMerger
from scripts.Validator.sequencevalidator import SequenceValidator
from scripts.MetaDataTable.metadatatable import MetadataTable


class MGExtract(SequenceValidator):
	"""
	Extracting Marker genes from fasta formatted files.
	"""

	_label = "MGExtract"

	_suffixes = {
		"5S": "5S_rRNA",
		"16S": "16S_rRNA",
		"23S": "23S_rRNA",
		"8S": "8S_rRNA",
		"18S": "18S_rRNA",
		"28S": "28S_rRNA"}

	_min_lengths = {
		"5S": 100,
		"16S": 900,
		"23S": 1000,
		"8S": 100,
		"18S": 900,
		"28S": 1000}

	def __init__(
		self, mg_analyse_executable, file_path_query_genome_file_paths, file_path_reference_genome_file_paths,
		file_path_name_reference_marker_genes, config_path, file_path_map_reference_genome_id_to_tax_id=None, max_processors=1, temp_directory=None,
		separator="\t", logfile=None, verbose=False, debug=False):
		"""
		Constructor

		@param mg_analyse_executable: File path to modified tool of Ivan
		@type mg_analyse_executable: str | unicode
		@param file_path_query_genome_file_paths: File path to file with the location of genomes to be classified
		@type file_path_query_genome_file_paths: str | unicode
		@param file_path_reference_genome_file_paths: File path to file with the location of reference genomes
		@type file_path_reference_genome_file_paths: str | unicode
		@param file_path_name_reference_marker_genes: File path to fasta file with list of marker gene sequences
		@type file_path_name_reference_marker_genes: str | unicode
		@param config_path: File path to configuration file
		@type config_path: str | unicode
		@param file_path_map_reference_genome_id_to_tax_id: Mapping of Reference genome_id to their taxonomic assignment
		@type file_path_map_reference_genome_id_to_tax_id: str | unicode
		@param max_processors: Amount of available processors
		@type max_processors: int | long
		@param temp_directory: File path to temporary storage
		@type temp_directory: str | unicode
		@param separator: Separator of metadata files
		@type separator: str | unicode
		"""
		super(MGExtract, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		assert file_path_map_reference_genome_id_to_tax_id is None or self.validate_file(file_path_map_reference_genome_id_to_tax_id)
		assert self.validate_file(file_path_query_genome_file_paths)
		assert file_path_reference_genome_file_paths is None or self.validate_file(file_path_reference_genome_file_paths)
		assert file_path_name_reference_marker_genes is None or self.validate_file(file_path_name_reference_marker_genes)
		assert self.validate_file(config_path)
		assert self.validate_file(mg_analyse_executable, executable=True)
		assert self.validate_number(max_processors, minimum=1)
		assert self.validate_dir(temp_directory)
		self._temp_directory = temp_directory
		self._mg_analyse_executable = mg_analyse_executable
		self._file_path_query_genome_file_paths = file_path_query_genome_file_paths
		self._file_path_reference_genome_file_paths = file_path_reference_genome_file_paths
		self._file_path_reference_marker_genes = file_path_name_reference_marker_genes
		self._config_path = config_path
		self._max_processors = max_processors
		self._debug = debug
		self._working_dirs = {}
		self._genome_id_to_tax_id = None
		self._separator = separator
		if file_path_map_reference_genome_id_to_tax_id is None:
			return
		meta_data_table = MetadataTable(
			separator=self._separator,
			logfile=logfile,
			verbose=verbose)
		meta_data_table.read(file_path_map_reference_genome_id_to_tax_id, column_names=False)
		self._genome_id_to_tax_id = meta_data_table.get_map(0, 1)
		del meta_data_table

	def gather_markergenes(self, hmmer, mg_type, file_path_output, file_path_map_uid_sid):
		"""
		Find and extract marker genes from genomes

		@param hmmer: hmmer2 or hmmer3
		@type hmmer: int | long
		@param mg_type: '16S', '5S' or '23S' etc
		@type mg_type: str | unicode
		@param file_path_output: Output for list of extracted marker genes sequences in fasta format
		@type file_path_output: str | unicode

		@rtype: None
		"""
		assert isinstance(hmmer, (int, long))
		assert isinstance(file_path_output, basestring)
		assert self.validate_number(hmmer, minimum=2, maximum=3)
		assert self.validate_dir(file_path_output, only_parent=True)
		assert mg_type in self._suffixes, "Marker gene '{}' is not supported."

		self._logger.info("Searching and extracting marker genes")
		start = time.time()
		query_genome_file_paths = self._get_genome_id_to_path_map(self._file_path_query_genome_file_paths)
		if self._file_path_reference_genome_file_paths is not None and self._file_path_reference_marker_genes is None:
			reference_genome_file_paths = self._get_genome_id_to_path_map(self._file_path_reference_genome_file_paths)
			query_genome_file_paths.update(reference_genome_file_paths)
		elif self._file_path_reference_genome_file_paths is not None and self._file_path_reference_marker_genes is not None:
			self._logger.warning("Ignoring reference genome file paths and using previous reference marker genes!")

		cmd_list = self._get_cmd_list(hmmer=hmmer, dict_of_fasta=query_genome_file_paths)
		list_of_tasks = []
		for cmd in cmd_list:
			list_of_tasks.append(parallel.TaskCmd(cmd))

		fail_list = parallel.runCmdParallel(list_of_tasks, self._max_processors)
		if fail_list is not None:
			for message in parallel.reportFailedCmd(fail_list):
				self._logger.error(message)
			msg = "Extracting marker genes failed."
			self._logger.error(msg)
			raise OSError(msg)

		tmp_out_file_path = tempfile.mktemp(suffix="_accepted", dir=self._temp_directory)
		tmp_out_file_bin_path = tempfile.mktemp(suffix="_rejected", dir=self._temp_directory)

		self._merge_marker_genes_files(
			query_genome_file_paths, tmp_out_file_path,
			file_path_out_bin=tmp_out_file_bin_path, file_path_map_uid_sid=file_path_map_uid_sid, mg_type=mg_type)
		if os.path.exists(tmp_out_file_path):
			shutil.copy2(tmp_out_file_path, file_path_output)
		else:
			self._logger.warning("No valid maker gene found!")
		if os.path.exists(tmp_out_file_bin_path):
			shutil.copy2(tmp_out_file_bin_path, file_path_output + ".rejected.fna")

		if self._file_path_reference_marker_genes is not None:
			# append reference genome marker genes
			shutil.copy(file_path_output, file_path_output + ".no_ref")
			with open(file_path_output, 'a') as write_handler, open(self._file_path_reference_marker_genes) as read_handler:
				write_handler.writelines(read_handler)

		end = time.time()
		self._logger.info("Extracting marker genes finished ({}s)".format(round(end - start, 1)))

		if not self._debug:
			for directory in self._working_dirs.values():
				shutil.rmtree(directory)
		else:
			for directory in self._working_dirs.values():
				self._logger.warning("Remove manually: '{}'".format(directory))

	def _get_cmd_list(self, hmmer, dict_of_fasta):
		"""
		Get list of command tasks

		@param hmmer: hmmer2 or hmmer3
		@type hmmer: int | long
		@param dict_of_fasta: Dict of file paths of genomes
		@type dict_of_fasta: dict[str|unicode, str|unicode]

		@return: List of command tasks
		@rtype: list[TaskCmd]
		"""
		assert isinstance(hmmer, (int, long))
		assert isinstance(dict_of_fasta, dict)
		assert self.validate_number(hmmer, minimum=2, maximum=3)
		arguments = [
			"-c '{}'".format(self._config_path),
			"-hmmer {}".format(hmmer),
			"-n",
		]
		cmd = "{exe} -i '{input_file}' -out '{dir_out}' {args}"
		if self._logfile:
			arguments.append("-log '{}'".format(self._logfile))
			cmd += " >> {}".format(self._logfile)

		cmd_list = []
		counter = 0
		for gid, file_path in dict_of_fasta.iteritems():
			self._working_dirs[gid] = tempfile.mkdtemp(dir=self._temp_directory, prefix="extract_{}".format(counter))
			cmd_list.append(cmd.format(
				exe=self._mg_analyse_executable, input_file=file_path, dir_out=self._working_dirs[gid], args=" ".join(arguments)))
			counter += 1
		self._logger.debug("\n".join(cmd_list))
		return cmd_list

	def _get_genome_id_to_path_map(self, file_path):
		"""
		Get a map of genome_id to genome path

		@param file_path: File path
		@type file_path: str | unicode

		@return: map of genome_id to genome path
		@rtype: dict[str|unicode, str|unicode]
		"""
		assert self.validate_file(file_path)

		data_table = MetadataTable(
			separator=self._separator,
			logfile=self._logfile,
			verbose=self._verbose)
		data_table.read(file_path, column_names=False)
		if data_table.get_number_of_rows() == 0:
			self._logger.warning("No data in file '{}'.".format(file_path))
			return {}
		dict_genome_id_to_path = data_table.get_map(0, 1)
		return dict_genome_id_to_path

	def _merge_marker_genes_files(self, dict_genome_id_to_path, file_path_out, file_path_out_bin, file_path_map_uid_sid, mg_type):
		"""
		Gather all marker gene sequences into single files

		@param dict_genome_id_to_path: Map of genome id to the real path of genomes
		@type dict_genome_id_to_path: dict[str|unicode, str|unicode]
		@param file_path_out: Fasta formatted output file path
		@type file_path_out: str|unicode
		@param file_path_out_bin: Fasta formatted output file path of rejected sequences
		@type file_path_out_bin: str|unicode
		@param mg_type: '16S', '5S' or '23S' etc
		@type mg_type: str | unicode

		@rtype: None
		"""
		assert mg_type in self._suffixes, "Marker gene '{}' is not supported."
		assert isinstance(dict_genome_id_to_path, dict)
		assert isinstance(file_path_out, basestring)
		assert isinstance(file_path_out_bin, basestring)
		suffix = self._suffixes[mg_type]
		min_length = self._min_lengths[mg_type]

		with open(file_path_out, "a") as stream_output,\
			open(file_path_out_bin, "a") as stream_output_bin,\
			open(file_path_map_uid_sid, "w") as stream_map_uid_sid:
			counter = 0
			sequence_merger = SequenceMerger(
				stream_output=stream_output,
				stream_output_bin=stream_output_bin,
				stream_map_uid_sid=stream_map_uid_sid,
				verbose=True)
			for genome_id, genome_path in dict_genome_id_to_path.iteritems():
				input_filename = os.path.basename(genome_path)
				input_filepath = "{prefix}.ids.{suffix}.fna".format(prefix=input_filename, suffix=suffix)
				input_file = os.path.join(self._working_dirs[genome_id], "working", input_filepath)
				# unique_id = genome_id
				unique_id = "GR_{}".format(counter)
				tax_id = ""
				if self._genome_id_to_tax_id is not None and genome_id in self._genome_id_to_tax_id:
					tax_id = self._genome_id_to_tax_id[genome_id]
				sequence_merger.merge(input_file, min_length=min_length, prefix_unique_id=unique_id, original_id=genome_id, tax_id=tax_id)
				counter += 1
				# concat_fasta_on_fasta.merge(
				# input_file, stream_output, min_length, unique_id=unique_id, out_bin_file_handle=stream_output_bin)
				# input_file_path, output_file_handle, min_length, unique_id=None, out_bin_file_handle=None, uid_sid_file_handle=None
