__author__ = 'hofmann'

import os
import sys
import argparse
import tempfile
from scripts.projectfilefolderhandle import ProjectFileFolderHandle
from scripts.configparserwrapper import ConfigParserWrapper
from scripts.Validator.validator import Validator


class ArgumentHandler(Validator):
	"""Reading pipeline configuration from file and from passed arguments"""
	_config_file_path = None
	_directory_pipeline = None
	_directory_temp = None

	# [main]
	_phase = 0
	_novelty_only = None
	_max_processors = 1

	# [MarkerGeneExtraction]
	_hmmer = None
	_file_path_reference_genome_locations = None
	_file_path_reference_markergene = None
	_file_path_query_genomes_location_file = None
	_file_path_map_reference_genome_id_to_tax_id = None
	_directory_output = None

	_binary_rnammer = None
	_hmmerBinDir = None  # 16S mg analysis
	_rnaHmmInstallDir = None  # 16S mg analysis
	_databaseFile = None  # 16S mg analysis

	# [MarkerGeneClustering]
	_cluster_method_choices = ['average', 'furthest', 'nearest']
	_binary_mothur = None
	metadata_table_in = None
	metadata_table_out = None
	cluster_method = None
	distance_cutoff = None
	silva_reference_directory = None
	precision = 1000

	# [MarkerGeneClassification]
	otu_distance = None
	classification_distance_minimum = None
	ncbi_reference_directory = None

	# [Binary]
	# binary_hmmer3 = None
	# binary_mummer = None

	# subfolder/files
	_silva_ref_files = ["mothur_ref_distances", "mothur_ref_names", "mothur_alignment_ref.fasta"]
	silva_ref_map_file = "map.tsv"
	_ncbi_ref_files = ["nodes.dmp", "merged.dmp", "names.dmp"]

	# meta table columns  'OTU', 'novelty_category'
	_separator = "\t"
	_column_name_genome_id = "genome_ID"
	_column_name_cutoff = "prediction_threshold"
	_column_name_otu_id = "OTU"
	_column_name_cluster_prediction = "NCBI_ID"
	_column_name_cluster_scientific_name = "SCIENTIFIC_NAME"
	_column_name_cluster_novelty = "novelty_category"
	_column_name_ani = "ANI"
	_column_name_ani_novelty = "ANI_NOVELTY_CATEGORY"
	_column_name_ani_compare = "ANI_TAXONOMIC_COMPARE"
	_column_name_ani_scientific_name = "ANI_SCIENTIFIC_NAME"

	def __init__(
		self, args=None, version="Prototype", separator="\t",
		column_name_genome_id="genome_ID", column_name_otu="OTU", column_name_novelty_category="novelty_category",
		column_name_ncbi="NCBI_ID"):

		self._separator = separator
		self._column_name_genome_id = column_name_genome_id
		self._column_name_otu_id = column_name_otu
		self._column_name_cluster_novelty = column_name_novelty_category
		self._column_name_ncbi = column_name_ncbi
		self._directory_pipeline = self._get_directory_pipeline()

		if not os.path.isabs(self._directory_pipeline):
			self._directory_pipeline = os.path.expanduser(self._directory_pipeline)
			self._directory_pipeline = os.path.realpath(self._directory_pipeline)

		self._valid_args = True

		# read parsed arguments
		self.parser = None
		options = self._get_parser_options(args, version)

		logfile = options.logfile
		if logfile is not None:
			logfile = self.get_full_path(logfile)
		super(ArgumentHandler, self).__init__(logfile=logfile)

		self._read_options(options)
		if not self._valid_args:
			return

		# set log level read from arguments
		self.set_log_level(verbose=self._verbose, debug=self._debug)

		# read config options
		self._read_config()
		if not self._valid_args:
			return

		# (sanity) check values
		self._check_values()

		tmp_dir = self._directory_temp
		directory_output = self._directory_output
		assert isinstance(tmp_dir, basestring)
		assert isinstance(directory_output, basestring)

		self._project_file_folder_handler = ProjectFileFolderHandle(
			tmp_dir=tmp_dir,
			output_dir=directory_output,
			time_stamp=None,
			logfile=self._logfile,
			verbose=self._verbose,
			debug=self._debug
		)

	def _get_directory_pipeline(self):
		"""
		Get pipeline location based on script location

		@return: Location of pipeline
		@rtype: str | unicode
		"""
		return self.get_full_path(os.path.dirname(os.path.realpath(sys.argv[0])))

	def to_file(self, file_path):
		file_directory = os.path.dirname(file_path)
		if not os.path.isdir(file_directory):
			self._logger.error("Directory does not exist: '{}'".format(file_directory))
			return
		with open(file_path, 'w') as file_handler:
			file_handler.write(self.to_string())

	def to_string(self):
		stages = ["Full", "MarkerGene Extraction", "MarkerGeneClustering", "MarkerGeneClustering", "ANI"]
		result_string = """Parameter:
		_Main_
		Config file:\t\t'{config}'
		Pipeline directory:\t'{pipe}'
		Output directory:\t'{out}'
		Stage:\t\t\t{stage}
		Processors:\t\t{pool}

		_MarkerGeneExtraction_
		Ref. genomes:\t\t'{ir}'
		Ref. 16S:\t\t'{irf}'
		New genomes:\t\t'{i}'

		_MarkerGeneClustering_
		Ref. SILVA:\t\t'{silva}
		Metadata Table in:\t'{im}
		Metadata Table out:\t'{om}
		Distance Cutoff:\t\t{th}

		_MarkerGeneClassification_
		Ref. NCBI:\t\t'{ncbi}'
		OTU dist.:\t\t{otu}
		Min. Clas. dist.:\t{mcd}

""".format(
			config=self._config_file_path,
			pipe=self._directory_pipeline,
			out=self._directory_output,
			stage=stages[self._phase],
			pool=self._max_processors,
			ir=self._file_path_reference_genome_locations,
			irf=self._file_path_reference_markergene,
			i=self._file_path_query_genomes_location_file,
			im=self.metadata_table_in,
			om=self.metadata_table_out,
			th=self.distance_cutoff,
			silva=self.silva_reference_directory,
			ncbi=self.ncbi_reference_directory,
			otu=self.otu_distance,
			mcd=self.classification_distance_minimum
		)
		return result_string

	def is_valid(self):
		return self._valid_args

	def _check_values(self):
		if self._novelty_only:
			if self._file_path_reference_genome_locations is None:
				self._logger.error("'-ir' Reference genome mapping file is required!")
				self._valid_args = False
				return
			elif not os.path.isfile(self._file_path_reference_genome_locations):
				self._logger.error("'-ir' File does not exist: '{}'".format(self._file_path_reference_genome_locations))
				self._valid_args = False
				return

			if self.metadata_table_out is None:
				self._logger.error("'-om' Metadata file is required!")
				self._valid_args = False
				return

			if self.metadata_table_in is None:
				self._logger.error("'-im' Metadata file is required!")
				self._valid_args = False
				return
			elif not os.path.isfile(self.metadata_table_in):
				self._logger.error("'-im' File does not exist: '{}'".format(self.metadata_table_in))
				self._valid_args = False
				return

			if self.ncbi_reference_directory is None:
				self._logger.error("'-ncbi' NCBI reference directory is required!")
				self._valid_args = False
				return

			self.ncbi_reference_directory = os.path.expanduser(self.ncbi_reference_directory)
			if not os.path.isabs(self.ncbi_reference_directory):
				self.ncbi_reference_directory = os.path.join(self._directory_pipeline, self.ncbi_reference_directory)

			if not os.path.isdir(self.ncbi_reference_directory):
				self._logger.error("'-ncbi' Directory does not exist: '{}'".format(self.ncbi_reference_directory))
				self._valid_args = False
				return

			for sub_directory in self._ncbi_ref_files:
				file_path = os.path.join(self.ncbi_reference_directory, sub_directory)
				if not os.path.isfile(file_path):
					self._logger.error("'-ncbi' File does not exist: '{}'".format(file_path))
					self._valid_args = False
					return
			return

		if self._directory_output is None and self.project_directory is None:
			self._logger.error("'-od' Output directory or '-o' project directory is required!")
			self._valid_args = False
			return

		if self._directory_output is not None and self.project_directory is None:
			directory_out = self._directory_output.rstrip('/')
			if not os.path.isdir(directory_out):
				directory_out = os.path.dirname(directory_out)

			if not os.path.isdir(directory_out):
				self._logger.error("'-od' Directory does not exist: '{}'".format(directory_out))
				self._valid_args = False
				return
			if self._phase >= 2:
				# search previous project
				list_of_items = os.listdir(self._directory_output)
				list_of_folder_directories = []
				for item in list_of_items:
					folder_item = os.path.join(self._directory_output, item)
					if os.path.isdir(folder_item):
						mg_fasta = self._project_file_folder_handler.get_file_path_mg_16s()
						if os.path.isfile(mg_fasta):
							list_of_folder_directories.append(folder_item)

				if len(list_of_folder_directories) == 0:
					self._logger.error("'-o' No valid Project directory found at '{}'".format(directory_out))
					self._valid_args = False
					return
				elif len(list_of_folder_directories) == 1:
					self.project_directory = list_of_folder_directories[0]
					self._logger.info("Selected project directory: '{}'".format(self.project_directory))
				else:
					self._logger.error("'-o' Several valid project directories found at '{}'. Please specify!".format(directory_out))
					self._valid_args = False
					return

		if not os.path.isdir(self.project_directory) and self._phase == 0:
			os.mkdir(self.project_directory)

		if self._directory_output is None:
			self._directory_output = os.path.dirname(self.project_directory.rstrip('/'))

		if not os.path.isabs(self.project_directory):
			self.project_directory = os.path.expanduser(self.project_directory)
			self.project_directory = os.path.realpath(self.project_directory)

		if not os.path.isabs(self._directory_output):
			self._directory_output = os.path.expanduser(self._directory_output)
			self._directory_output = os.path.realpath(self._directory_output)

		if not os.path.isdir(self._project_file_folder_handler.get_output_directory()):
			self._logger.error("'-o' Directory does not exist: '{}'".format(self._project_file_folder_handler.get_output_directory()))
			self._valid_args = False
			return

		self.silva_reference_directory = os.path.expanduser(self.silva_reference_directory)
		if self.silva_reference_directory is None:
			self._logger.error("'-silva' SILVA reference directory is required!")
			self._valid_args = False
			return
		if not os.path.isabs(self.silva_reference_directory):
			self.silva_reference_directory = os.path.realpath(self.silva_reference_directory)
		if not os.path.isdir(self.silva_reference_directory):
			self._logger.error("'-silva' Directory does not exist: '{}'".format(self.silva_reference_directory))
			self._valid_args = False
			return

		if self._directory_temp is None:
			self._directory_temp = tempfile.gettempdir()
		if not os.path.isdir(self._directory_temp):
			self._logger.error("'-temp' Directory does not exist: '{}'".format(self._directory_temp))
			self._valid_args = False
			return

		for sub_directory in self._silva_ref_files:
			file_path = os.path.join(self.silva_reference_directory, sub_directory)
			if not os.path.isfile(file_path):
				self._logger.error("'-silva' File does not exist: '{}'".format(file_path))
				self._valid_args = False
				return

		file_path = os.path.join(self.silva_reference_directory, self.silva_ref_map_file)
		if not os.path.isfile(file_path):
			self.silva_ref_map_file = None
			self._logger.warning("'-silva' A sequence mapping file for the silva references was not found. Assuming it is not required!")

		if self.ncbi_reference_directory is None:
			self._logger.error("'-ncbi' NCBI reference directory is required!")
			self._valid_args = False
			return

		self.ncbi_reference_directory = os.path.expanduser(self.ncbi_reference_directory)
		if not os.path.isabs(self.ncbi_reference_directory):
			self.ncbi_reference_directory = os.path.join(self._directory_pipeline, self.ncbi_reference_directory)

		if not os.path.isdir(self.ncbi_reference_directory):
			self._logger.error("'-ncbi' Directory does not exist: '{}'".format(self.ncbi_reference_directory))
			self._valid_args = False
			return

		for sub_directory in self._ncbi_ref_files:
			file_path = os.path.join(self.ncbi_reference_directory, sub_directory)
			if not os.path.isfile(file_path):
				self._logger.error("'-ncbi' File does not exist: '{}'".format(file_path))
				self._valid_args = False
				return

		if self._phase < 2:
			if self._databaseFile is None:
				self._logger.error("'databaseFile', bad argument!")
				self._valid_args = False
				return
			elif not os.path.isdir(self._databaseFile):
				self._logger.error("'databaseFile', bad path to directory! '{}'".format(self._databaseFile))
				self._valid_args = False
				return
			elif not os.path.isfile(os.path.join(self._databaseFile, "ncbitax_sqlite.db")):
				self._logger.error("'databaseFile', directory does not contain file: 'ncbitax_sqlite.db'! '{}'".format(self._databaseFile))
				self._valid_args = False
				return

			if not os.path.isdir(self._rnaHmmInstallDir):
				self._logger.error("'rnaHmmInstallDir', bad path to directory! '{}'".format(self._rnaHmmInstallDir))
				self._valid_args = False
				return

			if self._hmmer == 3 and self._hmmerBinDir is None:
				self._logger.error("'hmmerBinDir' is required for marker gene extraction, a path to the directory can be set in the config file!")
				self._valid_args = False
				return
			if self._hmmer == 2 and self._binary_rnammer is None:
				self._logger.error("'rnammer' is required for marker gene extraction, a path to the directory can be set in the config file!")
				self._valid_args = False
				return

			if self._hmmer == 3 and not os.path.isdir(self._hmmerBinDir):
				self._logger.error("'hmmerBinDir', bad path to directory! '{}'".format(self._hmmerBinDir))
				self._valid_args = False
				return
			if self._hmmer == 2 and not os.path.isfile(self._binary_rnammer):
				self._logger.error("'rnammer', bad path to binary! '{}'".format(self._binary_rnammer))
				self._valid_args = False
				return

			if self._hmmer == 3:
				executable = os.path.join(self._hmmerBinDir, "hmmsearch")
				rna_hmm_wrapper = os.path.join(self._rnaHmmInstallDir, "rna_hmm3.py")
			else:
				executable = self._binary_rnammer
				rna_hmm_wrapper = os.path.join(self._rnaHmmInstallDir, "rna_hmm2.py")

			if not os.access(executable, os.X_OK):
				self._logger.error("'hmmer{}', no permission to execute! '{}'".format(self._hmmer, executable))
				self._valid_args = False
				return
			if not os.access(rna_hmm_wrapper, os.X_OK):
				self._logger.error("'hmmer{}', no permission to execute! '{}'".format(self._hmmer, rna_hmm_wrapper))
				self._valid_args = False
				return

			if self._file_path_reference_genome_locations is None and self._file_path_reference_markergene is None:
				self._logger.error("'-ir' or '-irf' Reference genome maping file is required!")
				self._valid_args = False
				return
			elif self._file_path_reference_genome_locations is None and self._file_path_reference_markergene is None:
				self._logger.error("'-ir' or '-irf' Reference genome maping file is required!")
				self._valid_args = False
				return
			else:  # if not os.path.isfile(self.input_reference_file) and not os.path.isfile(self.input_reference_fna_file):
				file_path = self._file_path_reference_genome_locations or self._file_path_reference_markergene
				if not os.path.isfile(file_path):
					self._logger.error("'-ir','-irf' File does not exist: '{}'".format(file_path))
					self._valid_args = False
					return

		if self._phase <= 2:
			if self._binary_mothur is None:
				self._logger.error("'mothur' is required for clustering, a path to the binary can be set in in the config file!")
				self._valid_args = False
				return
			if not self.validate_file(self._binary_mothur, executable=True):
				self._valid_args = False
				return

		if self._phase > 2:
			cluster_file = self._project_file_folder_handler.get_file_path_cluster_mg_16s()
			if not os.path.isfile(cluster_file):
				self._logger.error("Mothur file with list of clusters not found at: '{}'".format(cluster_file))
				self._valid_args = False
				return

		if self._file_path_query_genomes_location_file is None:
			self._logger.error("'-i' Unidentified genome mapping file is required!")
			self._valid_args = False
			return
		elif not os.path.isfile(self._file_path_query_genomes_location_file):
			self._logger.error("'-i' File does not exist: '{}'".format(self._file_path_query_genomes_location_file))
			self._valid_args = False
			return

		if self.metadata_table_in is None:
			self._logger.error("'-im' Metadata file is required!")
			self._valid_args = False
			return
		elif not os.path.isfile(self.metadata_table_in):
			self._logger.error("'-im' File does not exist: '{}'".format(self.metadata_table_in))
			self._valid_args = False
			return

		if self._max_processors is None:
			self._logger.error("'-p' A number of processors is required!")
			self._valid_args = False
			return
		elif not self._max_processors > 0:
			self._logger.error("'-p' The number of processors must be a positive number: '{}'".format(self._max_processors))
			self._valid_args = False
			return

		if self.distance_cutoff is None:
			self._logger.error("'-th' A max distance threshold is required!")
			self._valid_args = False
			return
		elif not self.distance_cutoff > 0 or not self.distance_cutoff <= 1:
			self._logger.error("'-th' The max distance threshold must be between 0 and 1".format(self.distance_cutoff))
			self._valid_args = False
			return

		if self.otu_distance is None:
			self._logger.error("'-otu' A threshold is required!")
			self._valid_args = False
			return
		elif not self.otu_distance > 0 or not self.otu_distance <= 1:
			self._logger.error("'-otu' The number of processors must be a positive number: '{}'".format(self.otu_distance))
			self._valid_args = False
			return

		if self.classification_distance_minimum is None:
			self._logger.error("'-cth' A minimum classification distance threshold is required!")
			self._valid_args = False
			return
		elif not self.classification_distance_minimum > 0 or not self.classification_distance_minimum <= 1:
			self._logger.error("'-cth' The minimum classification distance threshold must be between 0 and 1: '{}'".format(self.classification_distance_minimum))
			self._valid_args = False
			return

		if self.cluster_method is None:
			self._logger.error("'-cm' A clustering method must be choosen: {}!".format(', '.join(self._cluster_method_choices)))
			self._valid_args = False
			return

		if self.cluster_method not in self._cluster_method_choices:
			self._logger.error("'-cm' A clustering method must be choosen: {}!".format(', '.join(self._cluster_method_choices)))
			self._valid_args = False
			return

		expected_output_size = self._expected_output_size_in_giga_byte()
		expected_tmp_size = expected_output_size
		# if self.multiple_samples:
		# 	expected_tmp_size /= self.number_of_samples
		directory_tmp = self._project_file_folder_handler.get_tmp_wd()
		directory_out = self._project_file_folder_handler.get_output_directory()
		free_tmp_space = self._free_space_in_giga_bytes(directory_tmp)
		free_out_space = self._free_space_in_giga_bytes(directory_out)
		message = None
		if expected_tmp_size > free_tmp_space:
			message = "WARNING: The directory '{dir}' has only {size:.4f} GigaByte of free space left. But output will require about {out_size:.4f} Gigabyte of space!\n".format(
				dir=directory_tmp,
				size=free_tmp_space,
				out_size=expected_tmp_size)
		elif expected_output_size > free_out_space:
			message = "WARNING: The directory '{dir}' has only {size:.4f} GigaByte of free space left. But output will require about {out_size:.4f} Gigabyte of space!\n".format(
				dir=directory_out,
				size=free_out_space,
				out_size=expected_output_size)
		elif expected_output_size > 100:
			message = "The output will require about {} GigaByte.".format(expected_output_size)

		if message is not None:
			user_input = raw_input(message + " Are you sure you want to continue? [y/n]\n>").lower()
			do_loop = True
			while do_loop:
				if user_input == 'n' or user_input == 'no':
					self._valid_args = False
					return
				if user_input == 'y' or user_input == 'yes':
					do_loop = False
					continue
				user_input = raw_input("Please type 'n' to abort, or 'y' to continue:\n>").lower()
		return

	# read the configuration file
	def _read_config(self):
		if not os.path.isfile(self._config_file_path):
			self._logger.error("'-c' File does not exist: '{}'".format(self._config_file_path))
			self._valid_args = False
			return

		self._config = ConfigParserWrapper(self._config_file_path, logfile=None, verbose=True)
		sections = ["Main", "MarkerGeneExtraction", "MarkerGeneClustering", "MarkerGeneClassification"]
		missing_section = self._config.validate_sections(sections)
		if missing_section:
			self._logger.error("Missing section '{}' in the configuration file.".format(missing_section))
			self._valid_args = False
			return
		# if self.pipeline_directory is None:
		# 	self.pipeline_directory = self._get_config_value("main", "pipeline_dir")

		if self._novelty_only is None:
			self._novelty_only = self._config.get_value("Main", "novelty_only", is_boolean=True, silent=False)

		if self._novelty_only:
			if self._file_path_reference_genome_locations is None:
				self._file_path_reference_genome_locations = self._config.get_value("MarkerGeneExtraction", "input_reference_file")

			if self.metadata_table_in is None:
				self.metadata_table_in = self._config.get_value("MarkerGeneClustering", "metadata_table_in")

			if self.metadata_table_out is None:
				self.metadata_table_out = self._config.get_value("MarkerGeneClustering", "metadata_table_out")

			if self.ncbi_reference_directory is None:
				self.ncbi_reference_directory = self._config.get_value("MarkerGeneClassification", "ncbi_reference_directory")

			return

		if self._directory_temp is None:
			self._directory_temp = self._config.get_value("Main", "temp_directory", silent=False)

		if self._directory_output is None:
			self._directory_output = self._config.get_value("Main", "output_directory", silent=False)

		if self._max_processors is None:
			self._max_processors = self._config.get_value("Main", "processors", is_digit=True)

		if self._binary_rnammer is None:
			self._binary_rnammer = self._config.get_value("MarkerGeneExtraction", "rnammer")

		if self._hmmerBinDir is None:
			self._hmmerBinDir = self._config.get_value("MarkerGeneExtraction", "hmmerBinDir")

		if self._databaseFile is None:
			self._databaseFile = self._config.get_value("MarkerGeneExtraction", "databaseFile")

		if self._rnaHmmInstallDir is None:
			self._rnaHmmInstallDir = self._config.get_value("MarkerGeneExtraction", "rnaHmmInstallDir")

		if self._file_path_reference_genome_locations is None:
			self._file_path_reference_genome_locations = self._config.get_value("MarkerGeneExtraction", "input_reference_file")

		if self._file_path_reference_markergene is None:
			self._file_path_reference_markergene = self._config.get_value("MarkerGeneExtraction", "input_reference_fna_file")

		if self._hmmer is None:
			self._hmmer = self._config.get_value("MarkerGeneExtraction", "hmmer", is_digit=True)

		if self._file_path_query_genomes_location_file is None:
			self._file_path_query_genomes_location_file = self._config.get_value("MarkerGeneExtraction", "input_genomes_file")

		if self._binary_mothur is None:
			self._binary_mothur = self._config.get_value("MarkerGeneClustering", "mothur")

		if self.metadata_table_in is None:
			self.metadata_table_in = self._config.get_value("MarkerGeneClustering", "metadata_table_in")

		if self.silva_reference_directory is None:
			self.silva_reference_directory = self._config.get_value("MarkerGeneClustering", "silva_reference_directory")

		if self.cluster_method is None:
			self.cluster_method = self._config.get_value("MarkerGeneClustering", "cluster_method")

		if self.ncbi_reference_directory is None:
			self.ncbi_reference_directory = self._config.get_value("MarkerGeneClassification", "ncbi_reference_directory")

		if self.distance_cutoff is None:
			self.distance_cutoff = self._config.get_value("MarkerGeneClustering", "max_threshold", is_digit=True)

		if self.otu_distance is None:
			self.otu_distance = self._config.get_value("MarkerGeneClustering", "otu_distance", is_digit=True)

		if self.classification_distance_minimum is None:
			self.classification_distance_minimum = self._config.get_value("MarkerGeneClustering", "classification_distance", is_digit=True)

	@staticmethod
	def _free_space_in_giga_bytes(directory="/tmp"):
		if not os.path.isdir(directory):
			return 0
		statvfs = os.statvfs(directory)
		free_space = statvfs.f_frsize * statvfs.f_bfree
		return free_space / float(1024 * 1024 * 1024)

	@staticmethod
	def _expected_output_size_in_giga_byte():
		expected_output_size = 0
		return expected_output_size

	def _read_options(self, options):
		config_file = options.config_file
		if config_file is not None:
			if not os.path.isabs(config_file):
				# config_file = os.path.join(self.pipeline_directory, self.folder_name_tools, config_file)
				config_file = os.path.realpath(config_file)
			if not os.path.isfile(config_file):
				self._logger.error("File does not exist: '{}'".format(options.config_file))
				self._valid_args = False
				return
		self._config_file_path = config_file
		self._verbose = options.verbose
		self.debug_mode = options.debug_mode
		self._logging = options.logging
		self._phase = options.stage
		self._max_processors = options.processors
		self._novelty_only = options.novelty_only

		self._file_path_reference_genome_locations = options.input_reference_file
		self._file_path_reference_markergene = options.input_reference_fna_file
		self._file_path_query_genomes_location_file = options.input_genomes
		self._directory_output = options.output_directory
		self.metadata_table_in = options.metadata_table_in
		self.cluster_method = options.cluster_method
		self.distance_cutoff = options.threshold
		self.otu_distance = options.otu_distance
		self.classification_distance_minimum = options.classification_distance

	def _get_parser_options(self, args=None, version="Prototype"):
		"""
		Parsing of passed arguments.

		@param args: Passed arguemnts

		@return: any
		"""
		description = """
	#######################################
	#    MetagenomeSimulationPipeline     #
	#    Version {}#
	#######################################

	Pipeline for the extraction of marker genes, clustering and taxonomic classification""".format(version.ljust(25))
		self.parser = argparse.ArgumentParser(
			usage="python %(prog)s configuration_file_path",
			version="MetagenomeSimulationPipeline {}".format(version),
			description=description,
			formatter_class=argparse.RawTextHelpFormatter)
		self.parser.add_argument("-verbose", "--verbose", action='store_true', default=False, help="display more information!")
		self.parser.add_argument("-debug", "--debug_mode", action='store_true', default=False, help="activate DEBUG modus. tmp folders will not be deleted!")
		self.parser.add_argument("-log", "--logging", action='store_true', default=False, help="pipeline output will written to a log file")
		# self.parser.add_argument("-qc", "--quality_check", action='store_true', default=False, help="")
		self.parser.add_argument(
			"-p", "--processors", default=None, type=int,
			help="number of processors to be used. >40 recommended.")
		self.parser.add_argument("-s", "--stage", default=0, type=int, choices=[0, 1, 2, 3, 4], help='''available options: 0-4:
0 -> Full run through,
1 -> Marker gene extraction,
2 -> Gene alignment and clustering,
3 -> Classification of Genomes and novelty prediction
4 -> Average Nucleotide Identity calculation
Default: 0
''')
		self.parser.add_argument("-n", "--novelty_only", action='store_true', default=None, help='''apply novelty categorisation only''')
		# group = parser.add_mutually_exclusive_group()
		group_input = self.parser.add_argument_group("input/output")
		group_input.add_argument(
			"-ir", "--input_reference_file", default=None, type=str,
			help="""path to a file containing list of reference genomes
Format: <genome_id>\\t<path>
No column names!""")
		group_input.add_argument(
			"-irf", "--input_reference_fna_file", default=None, type=str,
			help="path to a fasta file containing the 16S marker genes of the reference genomes")
		group_input.add_argument(
			"-i", "--input_genomes", default=None, type=str,
			help="""path to a file containing list of unidentified genomes
Format: <genome_id>\\t<path>
No column names!""")
		group_input.add_argument("-c", "--config_file", type=str, default=None, help="path to the configuration file of the pipeline")
		group_out = group_input.add_mutually_exclusive_group()
		group_out.add_argument(
			"-T", "--temp_directory", default=None, type=str,
			help="directory containing found marker genes and also a file in mothur format containing the clustering")
		group_out.add_argument(
			"-o", "--output_directory", default=None, type=str,
			help="directory containing found marker genes and also a file in mothur format containing the clustering")
		group_input.add_argument(
			"-im", "--metadata_table_in", default=None, type=str,
			help="path to file containing tab separated list of unidentified genomes")
		group_input.add_argument(
			"-om", "--metadata_table_out", default=None, type=str,
			help="path to file containing tab separated list of genomes and their file path")

		group_clustering = self.parser.add_argument_group("clustering")
		group_clustering.add_argument(
			"-cm", "--cluster_method", default=None,
			choices=self._cluster_method_choices, type=str,
			help="Algorithm used for clustering")
		group_clustering.add_argument(
			"-th", "--threshold", default=None, type=float,
			help="only distances up to the threshold will be calculated. Default: 0.05")
		group_clustering.add_argument(
			"-otu", "--otu_distance", default=None, type=float,
			help="genetic distances at which cluster will be used as otus. Default: 0.03")
		group_clustering.add_argument(
			"-cth", "--classification_distance", default=None, type=float,
			help="minimum distance for classification. Default: 0.02")

		if args is None:
			return self.parser.parse_args()
		else:
			return self.parser.parse_args(args)
