__author__ = 'Peter Hofmann'

import os
import shutil
import tempfile
import time
import scripts.MGExtract.concat_fasta_on_fasta as concat_fasta_on_fasta
import scripts.parallel as parallel
from scripts.Validator.sequencevalidator import SequenceValidator
from scripts.MetaDataTable.metadatatable import MetadataTable


class MGExtract(SequenceValidator):

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
		file_path_name_reference_marker_genes, config_path, file_path_iid_map=None, max_processors=1, temp_directory=None,
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
		@param file_path_iid_map:
		@type file_path_iid_map:
		@param max_processors: Amount of available processors
		@type max_processors: int | long
		@param temp_directory: File path to temporary storage
		@type temp_directory: str | unicode
		@param separator: Separator of metadata files
		@type separator: str | unicode
		"""
		super(MGExtract, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		assert file_path_iid_map is None or self.validate_file(file_path_iid_map)
		assert self.validate_file(file_path_query_genome_file_paths)
		assert self.validate_file(file_path_reference_genome_file_paths)
		assert self.validate_file(file_path_name_reference_marker_genes)
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
		self._working_dir = tempfile.mkdtemp(dir=self._temp_directory)
		self._iid_map = None
		self._separator = separator
		if file_path_iid_map is None:
			return
		meta_data_table = MetadataTable(
			separator=self._separator,
			logfile=logfile,
			verbose=verbose)
		meta_data_table.read(file_path_iid_map, column_names=False)
		self._iid_map = meta_data_table.get_map(0, 1)
		del meta_data_table

	def gather_markergenes(self, hmmer, mg_type, output_file):
		"""
		Find and extract marker genes from genomes

		@param hmmer: hmmer2 or hmmer3
		@type hmmer: int | long
		@param mg_type: '16S', '5S' or '23S' etc
		@type mg_type: str | unicode
		@param output_file: Output for list of extracted marker genes sequences in fasta format
		@type output_file: str | unicode

		@rtype: None
		"""
		assert isinstance(hmmer, (int, long))
		assert isinstance(output_file, basestring)
		assert self.validate_number(hmmer, minimum=2, maximum=3)
		assert self.validate_dir(output_file, only_parent=True)
		assert mg_type in self._suffixes, "Marker gene '{}' is not supported."

		self._logger.info("Searching and extracting marker genes")
		start = time.time()
		query_genome_file_paths = self._get_genome_id_to_path_map(self._file_path_query_genome_file_paths)
		if self._file_path_reference_genome_file_paths is not None and self._file_path_reference_marker_genes is None:
			reference_genome_file_paths = self._get_genome_id_to_path_map(self._file_path_reference_genome_file_paths)
			query_genome_file_paths.update(reference_genome_file_paths)
		elif self._file_path_reference_genome_file_paths is not None and self._file_path_reference_marker_genes is not None:
			self._logger.warning("Ignoring reference genome file paths and using previous reference marker genes!")

		# establish symbolic link to fasta files
		local_genome_file_paths = self._get_local_genome_paths(query_genome_file_paths)
		for genome_id, src in query_genome_file_paths.iteritems():
			dst = local_genome_file_paths[genome_id]
			if not self.validate_file(src):
				if not self._debug:
					shutil.rmtree(self._working_dir)
				else:
					self._logger.warning("Remove manually: '{}'".format(self._working_dir))
				return False
			os.symlink(src, dst)

		list_of_tasks = self._get_cmd_task_list(hmmer=hmmer, list_of_fasta=local_genome_file_paths.values())
		fail_list = parallel.runCmdParallel(list_of_tasks, self._max_processors)
		if fail_list is not None:
			for message in parallel.reportFailedCmd(fail_list):
				self._logger.error(message)
			msg = "Extracting marker genes failed."
			self._logger.error(msg)
			raise OSError(msg)

		tmp_out_file_path = tempfile.mktemp(suffix="_accepted", dir=self._working_dir)
		tmp_out_file_bin_path = tempfile.mktemp(suffix="_rejected", dir=self._working_dir)

		self._merge_marker_genes_files(local_genome_file_paths, tmp_out_file_path, out_bin_file_path=tmp_out_file_bin_path, mg_type=mg_type)
		if os.path.exists(tmp_out_file_path):
			shutil.copy2(tmp_out_file_path, output_file)
		else:
			self._logger.warning("No valid maker gene found!")
		if os.path.exists(tmp_out_file_bin_path):
			shutil.copy2(tmp_out_file_bin_path, output_file + ".rejected.fna")

		if self._file_path_reference_marker_genes is not None:
			# append reference genome marker genes
			shutil.copy(output_file, output_file + ".no_ref")
			with open(output_file, 'a') as write_handler, open(self._file_path_reference_marker_genes) as read_handler:
				write_handler.writelines(read_handler)

		end = time.time()
		self._logger.info("Extracting marker genes finished ({}s)".format(round(end - start, 1)))

		if not self._debug:
			shutil.rmtree(self._working_dir)
		else:
			self._logger.warning("Remove manually: '{}'".format(self._working_dir))

	def _get_cmd_task_list(self, hmmer, list_of_fasta):
		"""
		Get list of command tasks

		@param hmmer: hmmer2 or hmmer3
		@type hmmer: int | long
		@param list_of_fasta: List of file paths of genomes
		@type list_of_fasta: list[str|unicode]

		@return: List of command tasks
		@rtype: list[TaskCmd]
		"""
		assert isinstance(hmmer, (int, long))
		assert isinstance(list_of_fasta, list)
		assert self.validate_number(hmmer, minimum=2, maximum=3)
		out_dir = self._working_dir
		cmd = "{exe} -c '{config}' -nn -hmmer {hmmer} -i '{input_file}' -out '{out_dir}'"
		cmd_list = [cmd.format(exe=self._mg_analyse_executable, config=self._config_path, hmmer=hmmer, input_file=file_path, out_dir=out_dir) for file_path in list_of_fasta]
		return [parallel.TaskCmd(cmd, out_dir) for cmd in cmd_list]

	def _get_local_genome_paths(self, dict_genome_id_to_path):
		"""
		Create local paths for genomes to avoid several issues with mothur.

		@param dict_genome_id_to_path: Map of genome id to the real path of genomes
		@type dict_genome_id_to_path: dict[str|unicode, str|unicode]

		@return: Map of genome id to the local path of genomes
		@rtype: dict[str|unicode, str|unicode]
		"""
		assert isinstance(dict_genome_id_to_path, dict)

		system_link_directory = os.path.join(self._working_dir, "sym_links")
		if not os.path.exists(system_link_directory):
			os.makedirs(system_link_directory)
		dict_genome_id_to_local_path = {}
		for genome_id in dict_genome_id_to_path:
			basename = os.path.basename(dict_genome_id_to_path[genome_id])
			dict_genome_id_to_local_path[genome_id] = os.path.join(system_link_directory, basename)
		return dict_genome_id_to_local_path

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

	def _merge_marker_genes_files(self, dict_genome_id_to_path, out_file_path, out_bin_file_path, mg_type):
		"""
		Gather all marker gene sequences into single files

		@param dict_genome_id_to_path: Map of genome id to the real path of genomes
		@type dict_genome_id_to_path: dict[str|unicode, str|unicode]
		@param out_file_path: Fasta formatted output file path
		@type out_file_path: str|unicode
		@param out_bin_file_path: Fasta formatted output file path of rejected sequences
		@type out_bin_file_path: str|unicode
		@param mg_type: '16S', '5S' or '23S' etc
		@type mg_type: str | unicode

		@rtype: None
		"""
		assert mg_type in self._suffixes, "Marker gene '{}' is not supported."
		assert isinstance(dict_genome_id_to_path, dict)
		assert isinstance(out_file_path, basestring)
		assert isinstance(out_bin_file_path, basestring)
		suffix = self._suffixes[mg_type]
		min_length = self._min_lengths[mg_type]

		with open(out_file_path, "a") as output_file_handle, open(out_bin_file_path, "a") as out_bin_file_handle:
			for genome_id, genome_path in dict_genome_id_to_path.iteritems():
				input_filename = os.path.basename(genome_path)
				input_filepath = "{prefix}.ids.{suffix}.fna".format(prefix=input_filename, suffix=suffix)
				input_file = os.path.join(self._working_dir, "working", input_filepath)
				unique_id = genome_id
				if self._iid_map is not None:
					unique_id = self._iid_map[genome_id]
				concat_fasta_on_fasta.merge(input_file, output_file_handle, min_length, unique_id=unique_id, out_bin_file_handle=out_bin_file_handle)
				# input_file_path, output_file_handle, min_length, unique_id=None, out_bin_file_handle=None, uid_sid_file_handle=None
