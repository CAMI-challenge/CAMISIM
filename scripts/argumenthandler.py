__author__ = 'hofmann'
__version__ = '0.0.3'

import os
import sys
import argparse
import tempfile
import random
import numpy.random as np_random
import scripts
from scripts.configparserwrapper import ConfigParserWrapper
from scripts.Validator.validator import Validator, DefaultLogging
from scripts.ComunityDesign.communitydesign import Community
from scripts.projectfilefolderhandle import ProjectFileFolderHandle


class DefaultValues(DefaultLogging):
	"""
	Reading and writing config file

	@type _list_of_communities: list[Community]
	"""
	_seed = None

	_phase = 0
	_phase_validate_raw_genomes = False  # TODO: read from config
	_phase_design_community = False
	_phase_move_and_clean_genomes = False
	_phase_simulate_reads = False
	_phase_gsa = False
	_phase_pooled_gsa = False
	_phase_anonymize = False
	_phase_compress = False

	_separator = None
	_compresslevel = 0

	# ############
	# executables
	# ############
	_executable_art_illumina = None
	_executable_samtools = None

	# ############
	# reference directories
	# ############
	_directory_art_error_profiles = None
	_directory_ncbi_taxdump = None

	# ############
	# [main]
	# ############
	_tmp_dir = None
	_directory_output = None
	_directory_pipeline = None
	_max_processors = 1
	_dataset_id = ''

	# ############
	# [read_simulator]
	# ############
	_sample_size_in_base_pairs = None

	_read_simulator_type = None
	_error_profile = None
	_fragment_size_standard_deviation_in_bp = None
	_fragments_size_mean_in_bp = None

	# ############
	# [sampledesign]
	# ############
	_strain_simulation_template = None  # "tools/sgEvolver/simulation_dir"
	_number_of_samples = None
	_file_path_plasmid_sequence_names = None

	# ############
	# [comdesign]
	# ############
	_list_of_communities = []
	_input_list_of_file_paths_distributions = None

	def __init__(self, logfile=None, verbose=False, debug=False):
		super(DefaultValues, self).__init__(label="ArgumentHandler", logfile=logfile, verbose=verbose, debug=debug)
		self._validator = Validator(logfile=logfile, verbose=verbose, debug=debug)
		original_wd = os.getcwd()
		pipeline_dir = self._validator.get_full_path(os.path.dirname(scripts.__file__))
		os.chdir(pipeline_dir)
		self._DEFAULT_seed = random.randint(0, 2147483640)

		self._DEFAULT_phase = 0
		self._DEFAULT_phase_validate_raw_genomes = False  # TODO: read from config
		self._DEFAULT_phase_design_community = True
		self._DEFAULT_phase_move_and_clean_genomes = True
		self._DEFAULT_phase_simulate_reads = True
		self._DEFAULT_phase_gsa = True
		self._DEFAULT_phase_pooled_gsa = True
		self._DEFAULT_phase_anonymize = True
		self._DEFAULT_phase_compress = True

		self._DEFAULT_compresslevel = 5

		# ############
		# executables
		# ############
		self._DEFAULT_executable = 'art_illumina'
		self._DEFAULT_executable_samtools = 'samtools'

		# ############
		# reference directories
		# ############
		self._DEFAULT_directory_art_error_profiles = None
		self._DEFAULT_directory_ncbi_taxdump = None

		# ############
		# [main]
		# ############
		self._DEFAULT_tmp_dir = tempfile.gettempdir()
		# self._DEFAULT_directory_output = tempfile.mkdtemp(prefix="Output", dir=pipeline_dir)
		self._DEFAULT_directory_pipeline = pipeline_dir
		self._DEFAULT_max_processors = 1
		self._DEFAULT_dataset_id = 'default'

		# ############
		# [read_simulator]
		# ############
		self._DEFAULT_sample_size_in_base_pairs = 10 * 1000000000

		self._DEFAULT_read_simulator_type = 'art'
		self._DEFAULT_error_profile = 'hi150'
		self._DEFAULT_fragment_size_standard_deviation_in_bp = '270'
		self._DEFAULT_fragments_size_mean_in_bp = '27'

		# ############
		# [sampledesign]
		# ############
		self._DEFAULT_strain_simulation_template = os.path.join(
			pipeline_dir, 'scripts', 'StrainSimulationWrapper', 'sgEvolver', 'simulation_dir')
		self._DEFAULT_number_of_samples = 5
		self._DEFAULT_file_path_plasmid_sequence_names = None

		# ############
		# [comdesign]
		# ############
		# self._DEFAULT_number_of_communities = None
		self._DEFAULT_input_list_of_file_paths_distributions = None

		# self._DEFAULT_file_path_metadata_table =
		# self._DEFAULT_file_path_genome_locations =
		# self._DEFAULT_file_path_gff_locations =
		# self._DEFAULT_genomes_total =
		# self._DEFAULT_genomes_real =
		self._DEFAULT_limit_per_otu = 3
		self._DEFAULT_ratio = 1
		self._DEFAULT_mode = 'differential'
		self._DEFAULT_log_mu = 1
		self._DEFAULT_log_sigma = 2
		self._DEFAULT_gauss_mu = 1
		self._DEFAULT_gauss_sigma = 1
		self._DEFAULT_view = False

		os.chdir(original_wd)

	def _set_default_values(self):
		self._seed = self._seed or self._DEFAULT_seed

		self._phase = self._phase or self._DEFAULT_phase
		self._phase_validate_raw_genomes = self._phase_validate_raw_genomes or self._DEFAULT_phase_validate_raw_genomes
		self._phase_design_community = self._phase_design_community or self._DEFAULT_phase_design_community
		self._phase_move_and_clean_genomes = self._phase_move_and_clean_genomes or self._DEFAULT_phase_move_and_clean_genomes
		self._phase_simulate_reads = self._phase_simulate_reads or self._DEFAULT_phase_simulate_reads
		self._phase_gsa = self._phase_gsa or self._DEFAULT_phase_gsa
		self._phase_pooled_gsa = self._phase_pooled_gsa or self._DEFAULT_phase_pooled_gsa
		self._phase_anonymize = self._phase_anonymize or self._DEFAULT_phase_anonymize
		self._phase_compress = self._phase_compress or self._DEFAULT_phase_compress

		self._compresslevel = self._compresslevel or self._DEFAULT_compresslevel

		# ############
		# executables
		# ############
		self._executable_art_illumina = self._executable_art_illumina or self._DEFAULT_executable
		self._executable_samtools = self._executable_samtools or self._DEFAULT_executable_samtools

		# ############
		# reference directories
		# ############
		self._directory_art_error_profiles = self._directory_art_error_profiles or self._DEFAULT_directory_art_error_profiles
		self._directory_ncbi_taxdump = self._directory_ncbi_taxdump or self._DEFAULT_directory_ncbi_taxdump

		# ############
		# [main]
		# ############
		self._tmp_dir = self._tmp_dir or self._DEFAULT_tmp_dir
		self._directory_pipeline = self._directory_pipeline or self._DEFAULT_directory_pipeline
		if self._directory_output is None:
			self._directory_output = tempfile.mkdtemp(prefix="Output", dir=self._directory_pipeline)
		self._max_processors = self._max_processors or self._DEFAULT_max_processors
		self._dataset_id = self._dataset_id or self._DEFAULT_dataset_id

		# ############
		# [read_simulator]
		# ############
		self._sample_size_in_base_pairs = self._sample_size_in_base_pairs or self._DEFAULT_sample_size_in_base_pairs

		self._read_simulator_type = self._read_simulator_type or self._DEFAULT_read_simulator_type
		self._error_profile = self._error_profile or self._DEFAULT_error_profile
		self._fragment_size_standard_deviation_in_bp = self._fragment_size_standard_deviation_in_bp or self._DEFAULT_fragment_size_standard_deviation_in_bp
		self._fragments_size_mean_in_bp = self._fragments_size_mean_in_bp or self._DEFAULT_fragments_size_mean_in_bp

		# ############
		# [sampledesign]
		# ############
		self._strain_simulation_template = self._strain_simulation_template or self._DEFAULT_strain_simulation_template
		self._number_of_samples = self._number_of_samples or self._DEFAULT_number_of_samples
		self._file_path_plasmid_sequence_names = self._file_path_plasmid_sequence_names or self._DEFAULT_file_path_plasmid_sequence_names

		# ############
		# [comdesign]
		# ############
		# self._number_of_communities = None
		# self._input_list_of_file_paths_distributions = None
		for community in self._list_of_communities:
			community.limit_per_otu = community.limit_per_otu or self._DEFAULT_limit_per_otu
			community.ratio = community.ratio or self._DEFAULT_ratio
			community.mode = community.mode or self._DEFAULT_mode
			community.log_mu = community.log_mu or self._DEFAULT_log_mu
			community.log_sigma = community.log_sigma or self._DEFAULT_log_sigma
			community.gauss_mu = community.gauss_mu or self._DEFAULT_gauss_mu
			community.gauss_sigma = community.gauss_sigma or self._DEFAULT_gauss_sigma
			community.verbose = community.verbose or self._DEFAULT_view



class ConfigFileHandler(DefaultValues):
	"""
	Reading and writing config file

	@type _list_of_communities: list[Community]
	"""
	# internal variables not set in config
	_file_name_config = "config.cfg"
	_ncbi_ref_files = ["nodes.dmp", "merged.dmp", "names.dmp"]
	_base_pairs_multiplication_factor = float(1000000000)  # 10**9

	def __init__(self, logfile=None, verbose=False, debug=False):
		super(ConfigFileHandler, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		self._validator = Validator(logfile=logfile, verbose=verbose, debug=debug)

	def _read_config(self, file_path_config):
		"""
		Read parameter from configuration file.

		@rtype: bool
		"""
		# TODO: check that all keys options make sense
		self._config = ConfigParserWrapper(file_path_config)
		if not self._validator.validate_file(file_path_config, key="Configuration file"):
			self._valid_args = False
			return

		# ##########
		# [Main]
		# ##########

		section = None  # "Main"
		if self._phase is None:
			self._phase = self._config.get_value("phase", is_digit=True)

		if self._seed is None:
			self._seed = self._config.get_value("seed")

		if self._max_processors is None:
			self._max_processors = self._config.get_value("max_processors", is_digit=True)

		if self._dataset_id is None:
			self._dataset_id = self._config.get_value("dataset_id")

		if self._directory_output is None:
			self._directory_output = self._config.get_value("output_directory", is_path=True)

		if self._tmp_dir is None:
			config_value = self._config.get_value("temp_directory", is_path=True)
			if config_value is not None:
				assert self._validator.validate_dir(config_value)
				self._tmp_dir = config_value

		self._phase_gsa = self._config.get_value("gsa", is_boolean=True)
		self._phase_pooled_gsa = self._config.get_value("pooled_gsa", is_boolean=True)

		self._compresslevel = self._config.get_value("compress", is_digit=True)

		self._phase_anonymize = self._config.get_value("anonymous", is_boolean=True)

		# ##########
		# [ReadSimulator]
		# ##########

		section = None  # "ReadSimulator"
		if self._sample_size_in_base_pairs is None:
			config_value = self._config.get_value("size", is_digit=True)
			if config_value is not None:
				self._sample_size_in_base_pairs = long(config_value * self._base_pairs_multiplication_factor)

		if self._read_simulator_type is None:
			self._read_simulator_type = self._config.get_value("type")

		if self._executable_samtools is None:
			self._executable_samtools = self._config.get_value("samtools", is_path=True)

		if self._executable_art_illumina is None:
			self._executable_art_illumina = self._config.get_value("art_illumina", silent=True, is_path=True)

		if self._directory_art_error_profiles is None:
			self._directory_art_error_profiles = self._config.get_value("art_error_profiles", silent=True, is_path=True)

		if self._error_profile is None:
			self._error_profile = self._config.get_value("profile")

		if self._fragment_size_standard_deviation_in_bp is None:
			self._fragment_size_standard_deviation_in_bp = self._config.get_value(
				"fragment_size_standard_deviation", is_digit=True)

		if self._fragments_size_mean_in_bp is None:
			self._fragments_size_mean_in_bp = self._config.get_value("fragments_size_mean", is_digit=True)

		# ##########
		# [CommunityDesign]
		# ##########

		if self._input_list_of_file_paths_distributions is None:
			input_list_of_file_paths_distributions = self._config.get_value("distribution_file_paths", is_path=True, silent=True)
			if input_list_of_file_paths_distributions is not None:
				self._input_list_of_file_paths_distributions = input_list_of_file_paths_distributions.split(',')

		section = None  # "CommunityDesign"
		if self._directory_ncbi_taxdump is None:
			self._directory_ncbi_taxdump = self._config.get_value("ncbi_taxdump", is_path=True)

		if self._strain_simulation_template is None:
			self._strain_simulation_template = self._config.get_value(
				"strain_simulation_template", silent=True, is_path=True)

		if self._number_of_samples is None:
			self._number_of_samples = self._config.get_value("number_of_samples", is_digit=True)

		# if self._number_of_communities is None:
		# 	self._number_of_communities = self._config.get_value('number_of_communities', is_digit=True)
		#
		# if self._number_of_communities is None:
		# 	self._logger.error("Bad number of communities!")
		# 	self._valid_arguments = False
		# 	return

		community_sections = set()
		community_key_options = {
			"genomes_total", 'genomes_real', 'max_strains_per_otu', 'ratio',
			'log_mu', 'log_sigma', 'gauss_mu', 'gauss_sigma'}
		for key_options in community_key_options:
			community_sections = community_sections.union(self._config.search_sections_of(key_options))

		self._list_of_communities = []
		is_valid = True
		for community_section in community_sections:
			file_path_metadata_table = self._config.get_value('metadata', community_section, is_path=True)
			file_path_genome_locations = self._config.get_value('id_to_genome_file', community_section, is_path=True)
			file_path_gff_locations = self._config.get_value('id_to_gff_file', community_section, is_path=True, silent=True)
			mode = self._config.get_value('mode', community_section)
			if not isinstance(file_path_metadata_table, basestring):
				is_valid = False
			if not isinstance(file_path_genome_locations, basestring):
				is_valid = False
			# if not isinstance(file_path_gff_locations, basestring):
			# 	is_valid = False
			if not isinstance(mode, basestring):
				is_valid = False

			if not is_valid:
				continue
			assert isinstance(file_path_metadata_table, basestring)
			assert isinstance(file_path_genome_locations, basestring)
			assert file_path_gff_locations is None or isinstance(file_path_gff_locations, basestring)
			assert isinstance(mode, basestring)
			new_community = Community(
				identifier=community_section,
				genomes_total=self._config.get_value('genomes_total', community_section, is_digit=True),
				genomes_real=self._config.get_value('genomes_real', community_section, is_digit=True),
				limit_per_otu=self._config.get_value('max_strains_per_otu', community_section, is_digit=True),
				file_path_metadata_table=file_path_metadata_table,
				file_path_genome_locations=file_path_genome_locations,
				file_path_gff_locations=file_path_gff_locations,
				ratio=self._config.get_value('ratio', community_section, is_digit=True),
				mode=mode,
				log_mu=self._config.get_value('log_mu', community_section, is_digit=True),
				log_sigma=self._config.get_value('log_sigma', community_section, is_digit=True),
				gauss_mu=self._config.get_value('gauss_mu', community_section, is_digit=True),
				gauss_sigma=self._config.get_value('gauss_sigma', community_section, is_digit=True),
				verbose=self._config.get_value('view', community_section, is_boolean=True)
			)
			self._list_of_communities.append(new_community)
			self._number_of_communities = len(self._list_of_communities)
		return is_valid

	def _stream_main(self, output_stream=sys.stdout):
		"""

		@param output_stream:
		"""
		output_stream.write("[Main]\n")
		output_stream.write("seed={}\n".format(self._seed or ""))
		output_stream.write("phase={}\n".format(self._phase))
		output_stream.write("max_processors={}\n".format(self._max_processors))
		output_stream.write("dataset_id={}\n".format(self._dataset_id))
		output_stream.write("output_directory={}\n".format(self._directory_output or ""))
		output_stream.write("temp_directory={}\n".format(self._tmp_dir or ""))
		output_stream.write("gsa={}\n".format(self._phase_gsa))
		output_stream.write("pooled_gsa={}\n".format(self._phase_pooled_gsa))
		output_stream.write("anonymous={}\n".format(self._phase_anonymize))
		output_stream.write("compress={}\n".format(self._compresslevel))

	def _stream_read_simulator(self, output_stream=sys.stdout):
		"""

		@param output_stream:
		"""
		output_stream.write("[ReadSimulator]\n")
		output_stream.write("art_illumina={}\n".format(self._executable_art_illumina))
		output_stream.write("art_error_profiles={}\n".format(self._directory_art_error_profiles or ""))
		output_stream.write("samtools={}\n".format(self._executable_samtools))
		output_stream.write("profile={}\n".format(self._error_profile))
		output_stream.write("size={}\n".format(self._sample_size_in_base_pairs/self._base_pairs_multiplication_factor))
		output_stream.write("type={}\n".format(self._read_simulator_type))
		output_stream.write("fragments_size_mean={}\n".format(self._fragments_size_mean_in_bp))
		output_stream.write("fragment_size_standard_deviation={}\n".format(self._fragment_size_standard_deviation_in_bp))

	def _stream_community_design(self, output_stream=sys.stdout):
		"""

		@param output_stream:
		"""
		output_stream.write("[CommunityDesign]\n")
		output_stream.write("distribution_file_paths={}\n".format(self._input_list_of_file_paths_distributions or ""))
		output_stream.write("ncbi_taxdump={}\n".format(self._directory_ncbi_taxdump or ""))
		output_stream.write("strain_simulation_template={}\n".format(self._strain_simulation_template or ""))
		output_stream.write("number_of_samples={}\n".format(self._number_of_samples))
		# output_stream.write("number_of_communities={}\n".format(self._number_of_communities))

	def _stream_communities(self, output_stream=sys.stdout):
		"""

		@param output_stream:
		"""
		for community in self._list_of_communities:
			output_stream.write("[{}]\n".format(community.id))
			output_stream.write("metadata={}\n".format(community.file_path_metadata_table))
			output_stream.write("id_to_genome_file={}\n".format(community.file_path_genome_locations or ""))
			output_stream.write("id_to_gff_file={}\n".format(community.file_path_gff_locations or ""))
			output_stream.write("genomes_total={}\n".format(community.genomes_total))
			output_stream.write("genomes_real={}\n".format(community.genomes_real))
			output_stream.write("max_strains_per_otu={}\n".format(community.limit_per_otu))
			output_stream.write("ratio={}\n".format(community.ratio))
			output_stream.write("mode={}\n".format(community.mode))
			output_stream.write("log_mu={}\n".format(community.log_mu))
			output_stream.write("log_sigma={}\n".format(community.log_sigma))
			output_stream.write("gauss_mu={}\n".format(community.gauss_mu))
			output_stream.write("gauss_sigma={}\n".format(community.gauss_sigma))
			output_stream.write("view={}\n".format(community.verbose))
			output_stream.write("\n")

	def write_config(self, file_path):
		with open(file_path, 'w') as write_handler:
			self._stream_main(write_handler)
			write_handler.write("\n")
			self._stream_read_simulator(write_handler)
			write_handler.write("\n")
			self._stream_community_design(write_handler)
			write_handler.write("\n")
			self._stream_communities(write_handler)


class ArgumentHandler(ConfigFileHandler):
	"""Reading pipeline configuration from file and from passed arguments"""

	_label = "ArgumentHandler"

	_separator = None
	_file_path_config = None

	_column_name_genome_id = "genome_ID",
	_column_name_otu = "OTU",
	_column_name_novelty_category = "novelty_category",
	_column_name_ncbi = "NCBI_ID",
	_column_name_source = "source",

	def __init__(
		self, args=None, version="Prototype", separator="\t",
		column_name_genome_id="genome_ID", column_name_otu="OTU", column_name_novelty_category="novelty_category",
		column_name_ncbi="NCBI_ID", column_name_source="source"):
		"""
		Constructor

		@param args: sys.args or Nothing
		@type args: list
		@param separator: Expected separator in metadata file
		@type separator: str | unicode
		@param column_name_genome_id: Column name of genome id
		@type column_name_genome_id: str | unicode
		@param column_name_otu: Column name of otu identifiers
		@type column_name_otu: str | unicode
		@param column_name_novelty_category: Column name of novelty category
		@type column_name_novelty_category: str | unicode
		@param column_name_ncbi: Column name of ncbi ids
		@type column_name_ncbi: str | unicode
		@param column_name_source: Column name of genome sources
		@type column_name_source: str | unicode

		@rtype: None
		"""
		self._separator = separator
		self._column_name_genome_id = column_name_genome_id
		self._column_name_otu = column_name_otu
		self._column_name_novelty_category = column_name_novelty_category
		self._column_name_ncbi = column_name_ncbi
		self._column_name_source = column_name_source

		options = self._get_parser_options(args, version)
		logfile = options.logfile
		if logfile is not None:
			logfile = self._validator.get_full_path(logfile)
		super(ArgumentHandler, self).__init__(logfile=logfile)
		self._directory_pipeline = self._get_directory_pipeline()

		self._valid_arguments = True

		# read passed arguments
		self._input_list_of_file_paths_distributions = None
		self._read_options(options)
		if not self._valid_arguments:
			return

		# set log level read from arguments
		self.set_log_level(verbose=self._verbose, debug=self._debug)

		# read configuration files
		self._valid_arguments = self._read_config(self._file_path_config)
		if not self._valid_arguments:
			return

		# set missing values with default values
		self._set_default_values()

		# sanity check values
		self._check_values()
		if not self._valid_arguments:
			return
		# options = ArgumentHandler(args)

		# example: tmp_dir = "/tmp"
		tmp_dir = self._tmp_dir
		directory_output = self._directory_output
		assert isinstance(tmp_dir, basestring)
		assert isinstance(directory_output, basestring)

		if self._seed is not None:
			random.seed(self._seed)
			np_random.seed(abs(hash(self._seed)) % 4294967295)  # numpy accepts only 32 bit integers

		assert isinstance(self._directory_output, basestring)
		self._project_file_folder_handler = ProjectFileFolderHandle(
			tmp_dir=tmp_dir,
			output_dir=directory_output,
			time_stamp=None,
			logfile=self._logfile,
			verbose=self._verbose,
			debug=self._debug
		)
		self._project_file_folder_handler.make_directory_structure(self._number_of_samples)
		self.write_config(os.path.join(self._project_file_folder_handler.get_output_directory(), self._file_name_config))

	def _get_directory_pipeline(self):
		"""
		Get pipeline location based on script location

		@return: Location of pipeline
		@rtype: str | unicode
		"""
		return self._validator.get_full_path(os.path.dirname(os.path.realpath(sys.argv[0])))

	def to_file(self, file_path):
		"""
		Save parameter to a file

		@param file_path: A file path
		@type file_path: str | unicode

		@rtype: None
		"""
		assert self._validator.validate_file(file_path)
		with open(file_path, 'w') as file_handler:
			file_handler.write(self.to_string())

	def to_string(self):
		"""
		Return parameter as string

		@return: parameter as string
		@rtype: str
		"""
		# expected_output_size = self._expected_output_size_in_giga_byte()
		# if expected_output_size < 0.001:
		# 	expected_output_size = 0.001
		result_string = """
[Main]
# Starting point of the simulation
# 0: Do it all, 1: Only community design, 2: Start with read simulation
phase={phase}

# Maximum number of available processors
max_processors={pool}

# ID prefix used for anonymous sequences
dataset_id={dataset}

# Directory where the output will be stored. Will be created if it does not exist.
output_directory={out}

# Directory for temporary data (Example: /tmp)
temp_directory={temp}

# Make perfect assembly based on simulated reads of each sample
gsa={gsa}

# Make perfect assembly based on simulated reads of all samples
pooled_gsa={pgsa}

# Anonymize all sequences
anonymous={anonymous}

# Compress output data (0-9) with 9 being the strongest compression, but very slow.
compress={compress}

# An optional seed to get consistent results
seed={seed}

[ReadSimulator]
# Samtools (http://www.htslib.org/) takes care of sam/bam files. Version 1.0 or higher required!
# file path to executable
samtools={samtools}

# ART_Illumina (2008-2015) (Weichun Huang at whduke@gmail.com). Version 2.3.6 recommended!
# file path to executable
art_illumina={art_illumina}

# Directory containing error profiles for ART_Illumina
art_error_profiles={error_profiles}

# Supported profiles: "mi": EmpMiSeq250R, "hi": EmpHiSeq2kR, "hi150": HiSeq2500L150R
profile={error}

# Simulate samples of this size (giga base pairs)
size={gbps}

# Read simulator type (only ART is currently supported)
type={readsim}

# Mean size (bp) of fragment simulated by ART (read length depends on error profile)
fragments_size_mean={fmean}
# Standard deviation (bp) of fragments simulated by ART
fragment_size_standard_deviation={fsd}



[CommunityDesign]
# Directory with files of a ncbi taxdump
ncbi_taxdump={ncbi_taxdump}

# Directory of strain simulation template (optional, required if trains are to be simulated)
strain_simulation_template={template}

# The amount of samples to be simulated
number_of_samples={samples}


""".format(
			# config=self._file_path_config,
			# pipe=self._directory_pipeline,
			phase=self._phase,
			pool=self._max_processors,
			dataset=self._dataset_id,
			out=self._directory_output,
			temp=self._tmp_dir,
			gsa=str(self._phase_pooled_gsa),
			pgsa=str(self._phase_pooled_gsa),
			anonymous=str(self._phase_anonymize),
			compress=self._phase_compress,
			seed=self._seed,
			# bps=self._sample_size_in_base_pairs,
			# out_size=expected_output_size,
			samtools=self._executable_samtools,
			art_illumina=self._executable_art_illumina,
			error_profiles=self._directory_art_error_profiles,
			error=self._error_profile,
			gbps=float(self._sample_size_in_base_pairs)/self._base_pairs_multiplication_factor,
			readsim=self._read_simulator_type,
			fmean=self._fragments_size_mean_in_bp,
			fsd=self._fragment_size_standard_deviation_in_bp,
			# plasmid=self.plasmid_file
			ncbi_taxdump=self._directory_ncbi_taxdump,
			template=self._strain_simulation_template,
			samples=self._number_of_samples,
			)

		for community in self._list_of_communities:
			com_string = """
[community{i}]
# Metadata table, required tabseparated columns: genome_ID, OTU, NCBI_ID, novelty_category
metadata={metadata}

# File with genome file locations. Format: GENOME_ID \t FILEPATH \n
id_to_genome_file={id_to_genome_file}

# File with genome gen annotation file locations. Format: GENOME_ID \t FILEPATH \n
id_to_gff_file={id_to_gff_file}

# Total number of genomes to be used based on this community
genomes_total={genomes_total}

# Real number of genomes to be drawn from this community, rest will be simulated
genomes_real={genomes_real}

# For more diversity, strains are drawn from the same otu only up to a maximum.
# Maximum is exceeded if no other genomes available.
max_strains_per_otu={max_strains_per_otu}

# Base pair ratio between communities.
# If one has a ratio of 1 and the other a ratio of 2, the second will have twice the genome abundance
ratio={ratio}

# Simulated distribution
# Options: replicates, timeseries_normal, timeseries_lognormal, differential
mode={mode}

# mu of a log distribution
log_mu={log_mu}
# sigma of a log distribution
log_sigma={log_sigma}

# mu of a gauss distribution
gauss_mu={gauss_mu}
# sigma of a gauss distribution
gauss_sigma={gauss_sigma}

# View and confirm distribution (requires x-window)
view={view}


""".format(
				i=community.id,
				metadata=community.file_path_metadata_table,
				id_to_genome_file=community.file_path_genome_locations,
				id_to_gff_file=community.file_path_gff_locations,
				genomes_total=community.genomes_total,
				genomes_real=community.genomes_real,
				max_strains_per_otu=community.limit_per_otu,
				ratio=community.ratio,
				mode=community.mode,
				log_mu=community.log_mu,
				log_sigma=community.log_sigma,
				gauss_mu=community.gauss_mu,
				gauss_sigma=community.gauss_sigma,
				view=community.verbose
			)

			result_string += com_string
		return result_string

	def is_valid(self):
		"""
		Returns True if all parameter seem valid.

		@return: True if all parameter seem valid
		@rtype: bool
		"""
		return self._valid_arguments

	# ###################
	# Sanity check values
	# ###################

	def _check_common_values(self):
		"""
		Validate parameter that are always required.

		@rtype: None
		"""
		if self._number_of_samples is None:
			self._logger.error("'-ns' No number of samples given!")
			self._valid_arguments = False
		elif not self._validator.validate_number(self._number_of_samples, minimum=1, key='-ns'):
			self._valid_arguments = False

		if self._number_of_communities is None or self._number_of_communities == 0:
			self._logger.error("Bad 'number of communities': {}".format(self._number_of_communities))
			self._valid_arguments = False
		elif not self._validator.validate_number(self._number_of_communities, minimum=1):
			self._valid_arguments = False

		# not sure about that
		if self._file_path_plasmid_sequence_names is not None and not self._validator.validate_file(
			self._file_path_plasmid_sequence_names, key=''):
			self._valid_arguments = False

	def _check_community_design_values(self):
		"""
		Validate parameter that are required for community design.

		@rtype: None
		"""
		for index in range(self._number_of_communities):
			community = self._list_of_communities[index]
			if not community.has_valid_values():
				self._logger.error("[community{index}] Has an invalid value!".format(index=index))
				self._valid_arguments = False

		if self._directory_ncbi_taxdump is None:
			self._logger.error("NCBI taxdump directory is required!")
			self._valid_arguments = False
		elif not self._validator.validate_dir(self._directory_ncbi_taxdump, file_names=self._ncbi_ref_files, key='-ncbi'):
			self._valid_arguments = False
		else:
			self._directory_ncbi_taxdump = self._validator.get_full_path(self._directory_ncbi_taxdump)

		if self._strain_simulation_template is not None and self._validator.validate_dir(self._strain_simulation_template):
			self._strain_simulation_template = self._validator.get_full_path(self._strain_simulation_template)

		if self._executable_samtools is None:
			self._logger.error("Samtools executable is required!")
			self._valid_arguments = False
		elif not self._validator.validate_file(self._executable_samtools, executable=True):
			self._valid_arguments = False
		else:
			self._executable_samtools = self._validator.get_full_path(self._executable_samtools)

	def _check_read_simulation_values(self):
		"""
		Validate parameter that are required for simulating reads.

		@rtype: None
		"""
		if self._dataset_id is None:
			self._dataset_id = ''

		if self._sample_size_in_base_pairs is None:
			self._logger.error("'-bp' A size in giga basepairs must be given!")
			self._valid_arguments = False
		elif not self._validator.validate_number(self._sample_size_in_base_pairs, minimum=0, key='-bp', zero=False):
			self._valid_arguments = False

		if self._read_simulator_type is None:
			self._logger.error("'-rs' No read simulator declared!")
			self._valid_arguments = False
		elif self._read_simulator_type == 'art' or self._read_simulator_type == 'wgsim':
			if self._directory_art_error_profiles is None:
				self._logger.error("Art illumina error profile directory is required!")
				self._valid_arguments = False
			elif not self._validator.validate_dir(self._directory_art_error_profiles):
				self._valid_arguments = False
			else:
				self._directory_art_error_profiles = self._validator.get_full_path(self._directory_art_error_profiles)

			if self._executable_art_illumina is None:
				self._logger.error("Art illumina executable is required!")
				self._valid_arguments = False
			elif not self._validator.validate_file(self._executable_art_illumina, executable=True):
				self._valid_arguments = False
			else:
				self._executable_art_illumina = self._validator.get_full_path(self._executable_art_illumina)

			if self._directory_art_error_profiles is None:
				self._logger.error("Art illumina error profile directory is required!")
				self._valid_arguments = False
			elif not self._validator.validate_dir(self._directory_art_error_profiles):
				self._valid_arguments = False
			else:
				self._directory_art_error_profiles = self._validator.get_full_path(self._directory_art_error_profiles)

			if self._error_profile is None:
				self._logger.error("'-ep' An error profile for 'art' was not chosen!")
				self._valid_arguments = False

			if self._fragments_size_mean_in_bp is None:
				self._logger.error("'-fmean' For the simulation with 'art' a mean size of the fragments is required!")
				self._valid_arguments = False
			elif not self._validator.validate_number(self._fragments_size_mean_in_bp, minimum=1, key='-fmean'):
				self._valid_arguments = False

			if self._fragment_size_standard_deviation_in_bp is None:
				self._logger.error("'-fsd' For the simulation with 'art' a standard_deviation of the fragments size is required!")
				self._valid_arguments = False
			elif not self._validator.validate_number(self._fragment_size_standard_deviation_in_bp, minimum=1, key='-fsd'):
				self._logger.error(
					"'-fsd' The standard_deviation of the fragments size must be a positive number: '{}'".format(
						self._fragment_size_standard_deviation_in_bp))
				self._valid_arguments = False
		else:
			self._logger.error("Only art illumina is currently supported!")

		expected_output_size = self._expected_output_size_in_giga_byte()
		expected_tmp_size = expected_output_size / self._number_of_samples
		assert isinstance(self._directory_output, basestring)
		directory_out = self._directory_output
		directory_tmp = self._tmp_dir
		if not os.path.isdir(directory_out):
			directory_out = os.path.dirname(directory_out)

		user_input_required = False
		if not self._validator.validate_free_space(directory_tmp, required_space_in_gb=expected_tmp_size):
			user_input_required = True
		elif not self._validator.validate_free_space(directory_out, required_space_in_gb=expected_output_size):
			user_input_required = True
		elif expected_output_size > 100:
			# message = "The output will require approximately {} GigaByte.".format(expected_output_size)
			self._logger.warning("The output will require approximately {} GigaByte.".format(expected_output_size))
		if user_input_required:
			if not self._verbose:
				self._logger.error("Continuation only possible with enabled user input!>")
				self._valid_arguments = False
				return

			message = "Are you sure you want to continue? [y/n]"
			if not self.get_confirmation(message):
				self._valid_arguments = False
				return

		if self._phase_compress:
			if self._compresslevel is None:
				self._logger.error("No compression level (0 - 9) given.")
				self._valid_arguments = False
			elif not self._validator.validate_number(self._compresslevel, 0, 9):
				self._valid_arguments = False

	def _check_values(self):
		"""
		Validate parameter.

		@rtype: None
		"""
		if self._max_processors is None:
			self._max_processors = 1

		if self._tmp_dir is None:
			self._tmp_dir = tempfile.gettempdir()
		elif not self._validator.validate_dir(self._tmp_dir, key="temp directory"):
			self._valid_arguments = False
		else:
			self._tmp_dir = self._validator.get_full_path(self._tmp_dir)

		subfolders = ["scripts"]

		if self._compresslevel > 0:
			self._phase_compress = True

		if self._directory_pipeline is None:
			self._logger.error("Pipeline directory is required!")
			self._valid_arguments = False
			return
		elif not self._validator.validate_dir(
			self._directory_pipeline,
			sub_directories=subfolders, key="pipeline directory"):
			self._valid_arguments = False
			return

		if self._directory_output is None:
			# self._logger.error("'-o' Output directory is required!")
			self._logger.error("Output directory is required!")
			self._valid_arguments = False
			return

		self._directory_output = self._validator.get_full_path(self._directory_output)
		if not self._validator.validate_dir(self._directory_output, only_parent=True, key='-o'):
			self._valid_arguments = False
			return

		if self._phase is None:
			self._phase = 0

		if self._phase == 0:
			self._phase_validate_raw_genomes = True
			self._phase_design_community = True
			self._phase_move_and_clean_genomes = True
			self._phase_simulate_reads = True

		if self._phase == 1:
			self._phase_validate_raw_genomes = True
			self._phase_design_community = True
			self._phase_move_and_clean_genomes = True

		if self._phase == 2:
			self._phase_simulate_reads = True

		if self._phase == 2 and not self._validator.validate_dir(self._directory_output, key='-o'):
			self._valid_arguments = False
			return

		self._check_common_values()

		if self._phase_design_community:
			self._check_community_design_values()

		if self._phase_simulate_reads:
			self._check_read_simulation_values()

	def _expected_output_size_in_giga_byte(self):
		"""
		Get rough estimation of output size in giga byte.

		@return: Output size in giga byte.
		@rtype: float
		"""
		# very rough estimation
		expected_output_size = self._sample_size_in_base_pairs / float(1000000000) * 3 * 2
		expected_output_size *= self._number_of_samples
		return expected_output_size

	def _read_options(self, options):
		"""
		Read passed arguments.

		@rtype: None
		"""
		if not self._validator.validate_file(options.config_file, key='-c'):
			self._valid_arguments = False
			return
		self._file_path_config = self._validator.get_full_path(options.config_file)
		self._verbose = options.verbose
		self._debug = options.debug_mode
		self._phase = options.phase
		self._dataset_id = options.data_set_id
		self._max_processors = options.max_processors
		self._seed = options.seed
		# self._directory_output = options.output_directory
		# self._sample_size_in_base_pairs = options.sample_size_gbp
		# if self._sample_size_in_base_pairs is not None:
		# 	self._sample_size_in_base_pairs = long(options.sample_size_gbp * self._base_pairs_multiplication_factor)
		# self.read_simulator = options.read_simulator
		# self._error_profile = options.error_profile
		# self._fragment_size_standard_deviation_in_bp = options.fragment_size_standard_deviation
		# self._fragments_size_mean_in_bp = options.fragments_size_mean
		# self.plasmid_file = options.plasmid_file
		# self._number_of_samples = options.number_of_samples
		# self._phase_pooled_gsa = options.pooled_gsa

	@staticmethod
	def _get_parser_options(args=None, version="Prototype"):
		"""
		Parsing of passed arguments.

		@param args: Passed arguemnts

		@return: any
		"""
		parser = argparse.ArgumentParser(
			usage="python %(prog)s configuration_file_path",
			version="MetagenomeSimulationPipeline {}".format(version),
			description="""
	#######################################
	#    MetagenomeSimulationPipeline     #
	#    Version {}#
	#######################################

	Pipeline for the simulation of a metagenome""".format(version.ljust(25)),
			formatter_class=argparse.RawTextHelpFormatter)

		parser.add_argument(
			"-verbose", "--verbose",
			action='store_true',
			default=False,  # set None if read from config
			help="display more information")
		parser.add_argument(
			"-debug", "--debug_mode",
			action='store_true',
			default=False,  # set None if read from config
			help="more information, also temporary data will not be deleted")
		parser.add_argument(
			"-log", "--logfile",
			default=None,
			type=str,
			help="output will also be written to this log file")

		group_input = parser.add_argument_group('optional config arguments')
		group_input.add_argument(
			"-seed",
			default=None,
			type=str,
			help="seed for random number generators")
		group_input.add_argument(
			"-s", "--phase",
			default=None,
			type=int,
			choices=[0, 1, 2],
			help='''available options: 0,1,2. Default: 0
0 -> Full run,
1 -> Only Comunity creation,
2 -> Only Readsimulator
''')
		group_input.add_argument(
			"-id", "--data_set_id",
			default=None,
			help="id of the dataset, part of prefix of read/contig sequence ids")
		group_input.add_argument(
			"-p", "--max_processors",
			default=None,
			type=int,
			help="number of available processors")

		# ##########
		# i/o
		# ##########

		group_input = parser.add_argument_group('required')
		group_input.add_argument(
			"config_file",
			type=str,
			help="path to the configuration file")

		# group_input.add_argument(
		# 	"-o", "--output_directory",
		# 	type=str,
		# 	default=None,
		# 	help="Directory in which the simulated data will be put into")

		# ##########
		# read simulator
		# ##########

# 		group_read_simulator = parser.add_argument_group('read simulator')
# 		group_read_simulator.add_argument(
# 			"-rs", "--read_simulator",
# 			type=str,
# 			default=None,
# 			choices=['art', 'pirs', 'pacbio'],
# 			help='''choose read simulator:
# 	'pacbio',
# 	'art' (Illumina),
# 	'pirs' (Illumina).
# Currently only 'art' is supported.''')
# 		group_read_simulator.add_argument(
# 			"-ep", "--error_profile",
# 			type=str,
# 			default=None,
# 			choices=['mi', 'hi', 'hi150'],
# 			help='''Art Illumina error profiles:
# 	'mi': MiSeq 250bp,
# 	'hi': HiSeq 100bp,
# 	'hi150': HiSeq 150bp''')
# 		group_read_simulator.add_argument(
# 			"-bp", "--sample_size_gbp",
# 			type=float,
# 			default=None,
# 			help="The approximate total size the output will have in giga base pairs")
# 		group_read_simulator.add_argument(
# 			"-fmean", "--fragments_size_mean",
# 			type=int,
# 			default=None,
# 			help='''Mean size of the fragments in base pairs.''')
# 		group_read_simulator.add_argument(
# 			"-fsd", "--fragment_size_standard_deviation",
# 			type=int,
# 			default=None,
# 			help="Standard deviation from the mean size of fragments in base pairs.")
#
# 		# ##########
# 		# sample design
# 		# ##########
#
# 		group_community_design = parser.add_argument_group('sample design')
# 		group_community_design.add_argument(
# 			"-ns", "--number_of_samples",
# 			type=int,
# 			default=None,
# 			help='''Number of samples to be simulated''')
# 		group_community_design.add_argument(
# 			"-pgsa", "--pooled_gsa",
# 			default=None,
# 			action='store_true',
# 			help="Reads of all samples are pooled and assembled")
#
		if args is None:
			return parser.parse_args()
		else:
			return parser.parse_args(args)

	def get_confirmation(self, message):
		"""
		Get user confirmation of a question

		@param message: Yes, No question.
		@type message: str | unicode

		@return: Answer of Question
		@rtype: bool
		"""
		user_input = raw_input("{}\n>".format(message)).lower()
		while True:
			if self._validator.is_boolean_state(user_input):
				return self._validator.get_boolean_state(user_input)
			user_input = raw_input("Please type 'n' for no, or 'y' for yes:\n>").lower()
