__author__ = 'hofmann'

import os
import sys
import argparse
import tempfile
from scripts.projectfilefolderhandle_ga import ProjectFileFolderHandle
from scripts.configparserwrapper import ConfigParserWrapper
from scripts.Validator.sequencevalidator import SequenceValidator
from scripts.MGCluster.mgcluster import MGCluster
from scripts.MetaDataTable.metadatatable import MetadataTable


class ArgumentHandler(SequenceValidator):
	"""
	Reading pipeline configuration from file and from passed arguments
	"""
	_label = "ArgumentHandler"
	_file_path_config = None
	_directory_pipeline = None
	_directory_temp = None

	# [main]
	_phase = 0
	_max_processors = 1

	# [MarkerGeneExtraction]
	_hmmer = None
	_file_path_reference_genome_locations = None
	_file_path_reference_markergene = None
	_file_path_query_genomes_location_file = None
	_file_path_map_reference_genome_id_to_tax_id = None
	_directory_output = None

	_mg_analyse_executable = None
	_binary_rnammer = None
	_hmmerBinDir = None  # 16S mg analysis
	_rnaHmmInstallDir = None  # 16S mg analysis

	# [MarkerGeneClustering]
	_cluster_method_choices = MGCluster._cluster_method_choices
	_binary_mothur = None
	_metadata_table_in = None
	_cluster_method = None
	_distance_cutoff = None
	_silva_reference_directory = None
	_precision = 1000

	# [MarkerGeneAnnotation]
	_otu_distance = None
	_classification_distance_minimum = None
	_ncbi_reference_directory = None
	_ani_minimum_alignment = 0.8

	# subfolder/files
	_silva_ref_files = MGCluster._silva_ref_files
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
		"""
		Constructor

		@param args: Past arguments like sys.args
		@type args: list[]
		@param version: Version of main Program
		@type version: str|unicode
		@param separator: Expected separator for data tables
		@type separator: str|unicode
		@param column_name_genome_id: Column name of genome ids
		@type column_name_genome_id: str|unicode
		@param column_name_otu: Column name of otu ids
		@type column_name_otu: str|unicode
		@param column_name_novelty_category: Column name of genome novelty
		@type column_name_novelty_category: str|unicode
		@param column_name_ncbi: Column name of taxonomic classification
		@type column_name_ncbi: str|unicode
		"""
		assert args is None or isinstance(args, list)
		assert isinstance(version, str)
		assert isinstance(separator, str)
		assert isinstance(column_name_genome_id, str)
		assert isinstance(column_name_otu, str)
		assert isinstance(column_name_novelty_category, str)
		assert isinstance(column_name_ncbi, str)
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
		assert isinstance(tmp_dir, str)
		assert isinstance(directory_output, str)
		if not self.validate_dir(tmp_dir) or not self.validate_dir(directory_output, only_parent=True):
			self._valid_args = False
			return

		self._project_file_folder_handler = ProjectFileFolderHandle(
			tmp_dir=tmp_dir,
			output_dir=directory_output,
			time_stamp=None,
			logfile=self._logfile,
			verbose=self._verbose,
			debug=self._debug
		)

	def _get_mg_analyse_executable(self):
		"""
		Get the path of the marker gene analyse main script

		@return: File path
		@rtype: str|unicode
		"""
		return os.path.join(self._directory_pipeline, "rnahmm", "run.py")

	def _get_directory_pipeline(self):
		"""
		Get pipeline location based on script location

		@return: Location of pipeline
		@rtype: str | unicode
		"""
		return self.get_full_path(os.path.dirname(os.path.realpath(sys.argv[0])))

	def to_file(self, file_path):
		"""
		Write arguments as configuration file

		@param file_path:
		@type file_path: str | unicode

		@rtype: None
		"""
		assert self.validate_dir(file_path, only_parent=True)
		file_directory = os.path.dirname(file_path)
		if not os.path.isdir(file_directory):
			self._logger.error("Directory does not exist: '{}'".format(file_directory))
			return
		with open(file_path, 'w') as file_handler:
			file_handler.write(self.to_string())

	def to_string(self):
		"""
		Return arguments as string

		@return: arguments as string
		@rtype: str | unicode
		"""
		result_string = """Parameter:
		_Main_
		Config file:\t\t'{config}'
		Pipeline directory:\t'{pipe}'
		Output directory:\t'{out}'
		Phase:\t\t\t{stage}
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
			config=self._file_path_config,
			pipe=self._directory_pipeline,
			out=self._directory_output,
			stage=self._phase,
			pool=self._max_processors,
			ir=self._file_path_reference_genome_locations,
			irf=self._file_path_reference_markergene,
			i=self._file_path_query_genomes_location_file,
			im=self._metadata_table_in,
			om=self._project_file_folder_handler.get_file_path_meta_data_table(),
			th=self._distance_cutoff,
			silva=self._silva_reference_directory,
			ncbi=self._ncbi_reference_directory,
			otu=self._otu_distance,
			mcd=self._classification_distance_minimum
		)
		return result_string

	def _input_valid(self):
		"""
		Return True if input seems valid.

		@return: True if arguemnts valid
		@rtype: bool
		"""
		return self._valid_args

	def _validate_genome_ids(self):
		"""
		Validate genome ids

		@return:
		"""
		file_path_reference_genome_locations = self._file_path_reference_genome_locations
		file_path_query_genomes_location_file = self._file_path_query_genomes_location_file
		silva_reference_directory = self._silva_reference_directory
		assert isinstance(file_path_reference_genome_locations, str)
		assert isinstance(file_path_query_genomes_location_file, str)
		assert isinstance(silva_reference_directory, str)
		data_table_reference = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		data_table_reference.read(file_path_reference_genome_locations)
		reference_gids = data_table_reference.get_column(0)
		reference_gids_set = set(reference_gids)
		if not len(reference_gids) == len(reference_gids_set):
			self._valid_args = False
			self._logger.error("Reference genome ids are not unique")
			return

		data_table_query = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		data_table_query.read(file_path_query_genomes_location_file)
		query_gids = data_table_query.get_column(0)
		query_gids_set = set(query_gids)
		if not len(query_gids) == len(query_gids_set):
			self._valid_args = False
			self._logger.error("Query genome ids are not unique")
			return

		data_table_silva = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		file_path_silva_map = os.path.join(silva_reference_directory, MGCluster.get_file_name_of_map())
		data_table_silva.read(file_path_silva_map)
		silver_ids_set = set(data_table_silva.get_column(1))
		# silva ids are allowed to be not unique

		if not query_gids_set.isdisjoint(reference_gids_set):
			self._valid_args = False
			self._logger.error("Reference and query genomes ids must be unique!")
			return
		if not query_gids_set.isdisjoint(silver_ids_set):
			self._valid_args = False
			self._logger.error("Silva and query genomes ids must be unique!")
			return

	def _check_values(self):
		"""
		Validating input arguments

		@rtype: None
		"""
		if not self.validate_dir(self._directory_output, only_parent=True, key="Output directory"):
			self._valid_args = False
			return

		if not self.validate_dir(self._ncbi_reference_directory, file_names=self._ncbi_ref_files, key="NCBI reference directory"):
			self._valid_args = False
			return

		if self._directory_temp is None:
			self._directory_temp = tempfile.gettempdir()
		if not self.validate_dir(self._directory_temp, key="Temp directory"):
			self._valid_args = False
			return

		if self._phase < 2:
			if not self.validate_dir(self._rnaHmmInstallDir, key="rnaHmmInstallDir"):
				self._valid_args = False
				return

			if not self.validate_dir(self._rnaHmmInstallDir, file_names=["rna_hmm2.py", "rna_hmm3.py"]):
				self._valid_args = False
				return

			directory_rna_hmm = self._rnaHmmInstallDir
			assert isinstance(directory_rna_hmm, str)
			rna_hmm_wrapper = None
			if self._hmmer == 3:
				if not self.validate_dir(self._hmmerBinDir, file_names=["hmmsearch"]):
					self._valid_args = False
					return
				directory = self._hmmerBinDir
				assert isinstance(directory, str)
				executable = os.path.join(directory, "hmmsearch")
				if not self.validate_file(executable, executable=True):
					self._valid_args = False
					return
				rna_hmm_wrapper = os.path.join(directory_rna_hmm, "rna_hmm3.py")
			elif self._hmmer == 2:
				if not self.validate_file(self._binary_rnammer, executable=True, key="rnammer"):
					self._valid_args = False
					return
				rna_hmm_wrapper = os.path.join(directory_rna_hmm, "rna_hmm2.py")

			if not self.validate_file(rna_hmm_wrapper, executable=True, key="hmmer{}".format(self._hmmer)):
				self._valid_args = False
				return

			if self._file_path_reference_genome_locations is None and self._file_path_reference_markergene is None:
				self._logger.error("'-ir' or '-irf' Reference genome maping file is required!")
				self._valid_args = False
				return
			else:  # if not os.path.isfile(self.input_reference_file) and not os.path.isfile(self.input_reference_fna_file):
				file_path = self._file_path_reference_genome_locations or self._file_path_reference_markergene
				if not self.validate_file(file_path, key="reference genome"):
					self._valid_args = False
					return

		if self._phase == 0 or self._phase > 1:
			if not self.validate_file(self._metadata_table_in, key="Metadata file"):
				self._valid_args = False
				return

			if not self.validate_dir(self._silva_reference_directory, file_names=self._silva_ref_files, key="SILVA reference directory"):
				self._valid_args = False
				return

			if not self.validate_file(self._binary_mothur, executable=True, key="mothur"):
				self._valid_args = False
				return

			if self._distance_cutoff is None:
				self._logger.error("A max distance threshold is required!")
				self._valid_args = False
				return
			elif not self.validate_number(self._distance_cutoff, minimum=0, maximum=1, zero=False, key="Max distance threshold"):
				self._valid_args = False
				return

			if self._otu_distance is None:
				self._logger.error("A threshold is required for otus!")
				self._valid_args = False
				return
			elif not self.validate_number(self._otu_distance, minimum=0, maximum=1, zero=False, key="OTU distance threshold"):
				self._valid_args = False
				return

			if self._classification_distance_minimum is None:
				self._logger.error("A minimum classification distance threshold is required!")
				self._valid_args = False
				return
			elif not self.validate_number(self._classification_distance_minimum, minimum=0, maximum=1, zero=False, key="Minimum classification distance threshold"):
				self._valid_args = False
				return

			if self._cluster_method is None:
				self._logger.error("A clustering method must be chosen: {}!".format(', '.join(self._cluster_method_choices)))
				self._valid_args = False
				return

			if self._cluster_method not in self._cluster_method_choices:
				self._logger.error("A clustering method must be chosen: {}!".format(', '.join(self._cluster_method_choices)))
				self._valid_args = False
				return

		if not self.validate_file(self._file_path_query_genomes_location_file, key="Query genome locations"):
			self._valid_args = False
			return

		if self._max_processors is None:
			self._logger.error("A number of available processors is required!")
			self._valid_args = False
			return
		elif not self.validate_number(self._max_processors, minimum=1, key="Available processors"):
			self._valid_args = False
			return

		expected_output_size_gb = self._expected_output_size_in_giga_byte()
		expected_tmp_size = expected_output_size_gb
		if not self.validate_free_space(directory=self._directory_temp, required_space_in_gb=expected_tmp_size):
			self._valid_args = False
			return

		directory_output = self._directory_output
		assert isinstance(directory_output, str)
		if not os.path.exists(directory_output):
			directory_output = os.path.dirname(directory_output)
		if not self.validate_free_space(directory=directory_output, required_space_in_gb=expected_output_size_gb):
			self._valid_args = False
			return

		if self._file_path_nucmer:
			self.validate_file(self._file_path_nucmer, executable=True)

		self._validate_genome_ids()

	# read the configuration file
	def _read_config(self):
		"""
		Read arguments from a configuration file

		@rtype: None
		"""
		if not self.validate_file(self._file_path_config, key="Configuration file"):
			self._valid_args = False
			return

		self._config = ConfigParserWrapper(self._file_path_config, logfile=self._logfile, verbose=self._verbose)
		section = "Main"
		if self._phase is None:
			self._phase = self._config.get_value("phase", is_digit=True, silent=False)
		self._directory_temp = self._config.get_value("temp_directory", is_path=True, silent=False)
		self._directory_output = self._config.get_value("output_directory", is_path=True, silent=False)
		if self._max_processors is None:
			self._max_processors = self._config.get_value("max_processors", is_digit=True)
		self._validate_genomes = self._config.get_value("validate_genomes", is_boolean=True)

		section = "MarkerGeneExtraction"
		self._binary_rnammer = self._config.get_value("rnammer", is_path=True)
		self._hmmerBinDir = self._config.get_value("hmmerBinDir", is_path=True)
		self._rnaHmmInstallDir = self._config.get_value("rnaHmmInstallDir", is_path=True)
		self._file_path_reference_genome_locations = self._config.get_value("reference_genomes_file", is_path=True)
		self._file_path_map_reference_genome_id_to_tax_id = self._config.get_value("reference_genomes_map_file", is_path=True)
		self._file_path_reference_markergene = self._config.get_value("input_reference_fna_file", is_path=True)
		self._hmmer = self._config.get_value("hmmer", is_digit=True)
		self._file_path_query_genomes_location_file = self._config.get_value("input_genomes_file", is_path=True)

		section = "MarkerGeneClustering"
		self._binary_mothur = self._config.get_value("mothur", is_path=True)
		self._metadata_table_in = self._config.get_value("metadata_table_in", is_path=True)
		self._silva_reference_directory = self._config.get_value("silva_reference_directory", is_path=True)
		self._cluster_method = self._config.get_value("cluster_method")
		self._distance_cutoff = self._config.get_value("max_threshold", is_digit=True)
		self._otu_distance = self._config.get_value("otu_distance", is_digit=True)
		self._classification_distance_minimum = self._config.get_value("classification_distance", is_digit=True)

		section = "MarkerGeneAnnotation"
		self._ncbi_reference_directory = self._config.get_value("ncbi_reference_directory", is_path=True)
		self._file_path_nucmer = self._config.get_value("nucmer", is_path=True)
		self._annotate_classify = self._config.get_value("classify", is_boolean=True)
		self._annotate_novelty = self._config.get_value("novelty", is_boolean=True)
		self._annotate_otu = self._config.get_value("otu", is_boolean=True)
		self._annotate_ani = self._config.get_value("ani", is_boolean=True)

	@staticmethod
	def _expected_output_size_in_giga_byte():
		"""
		Get expected output size of the data in giga bytes
		@todo: Write a good predicting algorithm

		@return: Expected output size of the data in giga bytes
		@rtype: int|long
		"""
		expected_output_size = 0
		return expected_output_size

	def _read_options(self, options):
		"""
		Read option from parsed arguments with argparse

		@param options: parser.parse_args()
		@type options: Any

		@rtype: None
		"""
		config_file = options.config_file
		if config_file is not None:
			config_file = self.get_full_path(config_file)
		self._file_path_config = config_file
		self._verbose = options.verbose
		self._debug = options.debug_mode
		self._phase = options.phase
		self._max_processors = options.max_processors

	@staticmethod
	def _get_parser_options(args=None, version="Prototype"):
		"""
		Parsing of passed arguments.

		@param args: Passed arguemnts

		@return: any
		"""
		description = """
	#######################################
	#    GenomeAnnotationPipeline         #
	#    Version {}#
	#######################################

	Pipeline for the extraction of marker genes, clustering and taxonomic classification""".format(version.ljust(25))
		parser = argparse.ArgumentParser(
			usage="python %(prog)s configuration_file_path",
			version="MetagenomeSimulationPipeline TC {}".format(version),
			description=description,
			formatter_class=argparse.RawTextHelpFormatter)
		parser.add_argument(
			"-verbose", "--verbose",
			action='store_true',
			default=False,
			help="display more information!")
		parser.add_argument(
			"-debug", "--debug_mode",
			action='store_true',
			default=False,
			help="tmp folders will not be deleted!")
		parser.add_argument(
			"-log", "--logfile",
			type=str,
			default=None,
			help="pipeline output will written to this log file")

		group_input = parser.add_argument_group('optional config arguments')
		group_input.add_argument(
			"-p", "--max_processors",
			default=None,
			type=int,
			help="number of available processors")
		group_input.add_argument("-s", "--phase", default=None, type=int, choices=[0, 1, 2, 3], help='''
0 -> Full run (Default)
1 -> Marker gene extraction
2 -> Gene alignment and clustering
3 -> Annotation of Genomes
''')
		group_input = parser.add_argument_group('required')
		group_input.add_argument("config_file", type=str, default=None, help="path to the configuration file of the pipeline")

		if args is None:
			return parser.parse_args()
		else:
			return parser.parse_args(args)
