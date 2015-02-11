__author__ = 'hofmann'

import os
import shutil
import tempfile
import time
import source.concat_fasta_on_fasta as concat_fasta_on_fasta
import source.parallel as parallel
from source.argumenthandler import ArgumentHandler
from source.logger import Logger


def main(options):
	assert isinstance(options, ArgumentHandler)
	mg_extract = MGExtract(mg_analyse_executable=os.path.join(options.pipeline_directory, "tools", "16SDetector", "run.py"),
						   filename_query_genome_file_paths=options.input_genomes_file,
						   filename_reference_genome_file_paths=options.input_reference_file,
						   filename_reference_marker_genes=options.input_reference_fna_file,
						   config_path=options.config_file_path,
						   max_processors=options.processors,
						   temp_directory=options.temp_directory,
						   debug=options._debug_mode)

	return mg_extract.gather_markergenes(hmmer=options.hmmer,
										 mg_type="16S",
										 output_file=os.path.join(options.project_directory, options.file_mg_16s))


class MGExtract(object):
	def __init__(self, mg_analyse_executable, filename_query_genome_file_paths, filename_reference_genome_file_paths,
				 filename_reference_marker_genes, config_path, max_processors=1, temp_directory=None, debug=False, logger=None):
		self._logger = logger
		if not self._logger:
			self._logger = Logger("MGExtract")
		self._temp_directory = temp_directory
		self._mg_analyse_executable = mg_analyse_executable
		self._filename_query_genome_file_paths = filename_query_genome_file_paths
		self._filename_reference_genome_file_paths = filename_reference_genome_file_paths
		self._filename_reference_marker_genes = filename_reference_marker_genes
		self._config_path = config_path
		self._max_processors = max_processors
		self._debug = debug
		self._working_dir = tempfile.mkdtemp(dir=self._temp_directory)

	def gather_markergenes(self, hmmer, mg_type, output_file):
		self._logger.info("[MGExtract] Searching and extracting marker genes")
		start = time.time()
		success = True
		query_genome_file_paths = self._parse_genome_file_path_file(self._filename_query_genome_file_paths)
		if self._filename_reference_genome_file_paths is not None and self._filename_reference_marker_genes is None:
			reference_genome_file_paths = self._parse_genome_file_path_file(self._filename_reference_genome_file_paths)
			query_genome_file_paths.update(reference_genome_file_paths)
		elif self._filename_reference_genome_file_paths is not None and self._filename_reference_marker_genes is not None:
			self._logger.warning("[MGExtract] Ignoring reference genome file paths and using previous reference marker genes!")
		local_genome_file_paths = self._get_local_genome_paths(query_genome_file_paths)

		# establish symbolic link to fasta files
		for genome_id, src in query_genome_file_paths.iteritems():
			dst = local_genome_file_paths[genome_id]
			if not os.path.exists(src):
				self._logger.error("[MGExtract] File does not exist: '{}'".format(src))
				if not self._debug:
					shutil.rmtree(self._working_dir)
				else:
					self._logger.warning("[MGExtract] Remove manually: '{}'".format(self._working_dir))
				return False
			os.symlink(src, dst)

		cmd_task_list = self._create_cmd_task_list(hmmer=hmmer, list_of_fasta=local_genome_file_paths.values())
		fail_list = parallel.runCmdParallel(cmd_task_list, self._max_processors)
		if fail_list is not None:
			success = False

		if success:
			tmp_out_file_path = tempfile.mktemp(suffix="_accepted", dir=self._temp_directory)
			tmp_out_file_bin_path = tempfile.mktemp(suffix="_rejected", dir=self._temp_directory)

			self._merge_marker_genes_files(local_genome_file_paths, tmp_out_file_path, out_bin_file=tmp_out_file_bin_path, mg_type=mg_type)
			if os.path.exists(tmp_out_file_path):
				shutil.copy2(tmp_out_file_path, output_file)
			if os.path.exists(tmp_out_file_bin_path):
				shutil.copy2(tmp_out_file_bin_path, output_file+".rejected.fna")

			if self._filename_reference_marker_genes is not None:
				# append reference genome marker genes
				shutil.copy(output_file, output_file+".no_ref")
				with open(output_file, 'a') as write_handler, open(self._filename_reference_marker_genes) as read_handler:
					write_handler.writelines(read_handler)
		end = time.time()
		self._logger.info("[MGExtract] Done ({}s)".format(round(end - start), 1))

		if not self._debug:
			shutil.rmtree(self._working_dir)
		else:
			self._logger.warning("[MGExtract] Remove manually: '{}'".format(self._working_dir))
		return success

	def _create_cmd_task_list(self, hmmer, list_of_fasta):
		out_dir = self._working_dir
		# TODO: use logfile instead of /dev/null, idealy detect errors
		#cmd = "{exe} -c '{config}' -nn -hmmer {hmmer} -i '{input_file}' -out '{out_dir}' > /dev/null 2> /dev/null"
		cmd = "{exe} -c '{config}' -nn -hmmer {hmmer} -i '{input_file}' -out '{out_dir}'"
		cmd_list = [cmd.format(exe=self._mg_analyse_executable, config=self._config_path, hmmer=hmmer, input_file=file_path, out_dir=out_dir) for file_path in list_of_fasta]
		return [parallel.TaskCmd(cmd, out_dir) for cmd in cmd_list]

	def _get_local_genome_paths(self, dict_genome_id_to_path):
		system_link_directory = os.path.join(self._working_dir, "sym_links")
		if not os.path.exists(system_link_directory):
			os.makedirs(system_link_directory)
		dict_genome_id_to_local_path = {}
		for genome_id in dict_genome_id_to_path:
			basename = os.path.basename(dict_genome_id_to_path[genome_id])
			dict_genome_id_to_local_path[genome_id] = os.path.join(system_link_directory, basename)
		return dict_genome_id_to_local_path

	@staticmethod
	def _parse_genome_file_path_file(file_path):
		dict_genome_id_to_path = {}
		with open(file_path) as read_handler:
			for line in read_handler:
				if line.startswith('#'):
					continue
				line = line.strip()
				if len(line) == 0:
					continue
				genome_id, genome_path = line.split('\t')
				dict_genome_id_to_path[genome_id] = genome_path
		return dict_genome_id_to_path

	# gather all marker gene sequences into single files
	def _merge_marker_genes_files(self, dict_genome_id_to_path, out_file_path, out_bin_file, mg_type):
		suffixes = {"5S": "5S_rRNA",
					"16S": "16S_rRNA",
					"23S": "23S_rRNA",
					"8S": "8S_rRNA",
					"18S": "18S_rRNA",
					"28S": "28S_rRNA"}

		min_lengths = {"5S": 100,
					"16S": 900,
					"23S": 1000,
					"8S": 100,
					"18S": 900,
					"28S": 1000}
		suffix = suffixes[mg_type]
		min_length = min_lengths[mg_type]
		assert isinstance(dict_genome_id_to_path, dict)
		for genome_id, genome_path in dict_genome_id_to_path.iteritems():
			input_filename = os.path.basename(genome_path)
			input_filepath = "{prefix}.ids.{suffix}.fna".format(prefix=input_filename,
																suffix=suffix)
			input_file = os.path.join(self._working_dir, "working", input_filepath)
			concat_fasta_on_fasta.merge(input_file, out_file_path, min_length, unique_id=genome_id, out_bin_file=out_bin_file)
