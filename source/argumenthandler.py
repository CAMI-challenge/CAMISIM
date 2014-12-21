__author__ = 'hofmann'

import os
import argparse
import time
import datetime
from config import Config
from logger import Logger


class ArgumentHandler(object):
	"""Reading pipeline configuration from file and from passed arguments"""
	_stages = 0
	_logging = False
	_debug_mode = False
	_verbose = True
	config_file_path = None
	pipeline_directory = None

	#[main]
	stage = 0
	processors = 1

	#[MarkerGeneExtraction]
	hmmer = None
	input_reference_file = None
	input_reference_fna_file = None
	input_genomes_file = None
	output_directory = None
	project_directory = None

	#[MarkerGeneClustering]
	metadata_table_in = None
	metadata_table_out = None
	distance_cutoff = None
	silva_reference_directory = None

	#[MarkerGeneClassification]
	otu_distance = None
	classification_distance_minimum = None
	ncbi_reference_directory = None

	#[Binary]
	binary_mothur = None
	binary_hmmer3 = None
	binary_rnammer = None
	binary_mummer = None

	#subfolder_names
	# name of folder containing all tools
	folder_name_source = "source"
	folder_name_tools = "tools"
	# name of folder containing all log files
	folder_name_logfiles = "logfiles"

	#filenames
	_filename_config_default = "config.cfg"
	_suffix_16S = "16S_rRNA"
	filename_log = "pipeline.log"
	#file_mg_23s = "23S_rRNA.fna"
	file_mg_16s = "{}.fna".format(_suffix_16S)
	#file_mg_05s = "5S_rRNA.fna"
	file_cluster_mg_16s = "mothur_cluster_{}.list".format(_suffix_16S)

	#subfolder/files
	_silva_ref_files = ["mothur_ref_distances", "mothur_ref_names", "mothur_alignment_ref.fasta"]
	_ncbi_ref_files = ["nodes.dmp", "merged.dmp", "names.dmp"]

	#meta table columns  'OTU', 'novelty_category'
	column_name_unpublished_genomes_id = "genome_ID"
	column_name_cutoff = "prediction_threshold"
	column_name_otu_id = "OTU"
	column_name_cluster_prediction = "NCBI_ID"
	column_name_cluster_scientific_name = "SCIENTIFIC_NAME"
	column_name_cluster_novelty = "novelty_category"
	column_name_ani = "ANI"
	column_name_ani_novelty = "ANI_NOVELTY_CATEGORY"
	column_name_ani_compare = "ANI_TAXONOMIC_COMPARE"
	column_name_ani_scientific_name = "ANI_SCIENTIFIC_NAME"

	def __init__(self, args=None, pipeline_directory=None, stages=1, logger=None):
		self._stages = stages

		if pipeline_directory is None:
			ArgumentHandler.pipeline_directory = os.getcwd()
		else:
			ArgumentHandler.pipeline_directory = pipeline_directory

		if not os.path.isabs(ArgumentHandler.pipeline_directory):
			ArgumentHandler.pipeline_directory = os.path.expanduser(ArgumentHandler.pipeline_directory)
			ArgumentHandler.pipeline_directory = os.path.realpath(ArgumentHandler.pipeline_directory)

		self._logger = logger
		if self._logger is None:
			self._logger = Logger("ArgumentHandler")

		self._valid_args = True

		# read parsed arguments
		self.parser = None
		options = self._get_parser_options(args)
		self._read_options(options)
		if not self._valid_args:
			return

		# read config options
		self._read_config()
		if not self._valid_args:
			return

		# (sanity) check values
		self._check_values()

	def set_logger(self, logger=None):
		if logger is None:
			self._logger.error("Could not set new logger!")
			return
		self._logger = logger

	def to_file(self, file_path):
		file_directory = os.path.dirname(file_path)
		if not os.path.isdir(file_directory):
			self._logger.error("Directory does not exist: '{}'".format(file_directory))
			return
		with open(file_path, 'w') as file_handler:
			file_handler.write(self.to_string())

	@staticmethod
	def to_string():
		stages = ["Full", "MarkerGene Extraction", "MarkerGeneClustering", "MarkerGeneClustering", "ANI"]
		result_string = """Parameter:
		_Main_
		Config file:\t\t'{config}'
		Pipeline directory:\t'{pipe}'
		Output directory:\t'{out}'
		Project Folder:\t\t'{project}'
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

""".format(	config=ArgumentHandler.config_file_path,
			pipe=ArgumentHandler.pipeline_directory,
			out=ArgumentHandler.output_directory,
			stage=stages[ArgumentHandler.stage],
			pool=ArgumentHandler.processors,
			ir=ArgumentHandler.input_reference_file,
			irf=ArgumentHandler.input_reference_fna_file,
			i=ArgumentHandler.input_genomes_file,
			project=ArgumentHandler.project_directory,
			im=ArgumentHandler.metadata_table_in,
			om=ArgumentHandler.metadata_table_out,
			th=ArgumentHandler.distance_cutoff,
			silva=ArgumentHandler.silva_reference_directory,
			ncbi=ArgumentHandler.ncbi_reference_directory,
			otu=ArgumentHandler.otu_distance,
			mcd=ArgumentHandler.classification_distance_minimum,
			)
		return result_string

	def is_valid(self):
		return self._valid_args

	def _check_values(self):
		if ArgumentHandler.output_directory is None and ArgumentHandler.project_directory is None:
			self._logger.error("'-od' Output directory or '-o' project directory is required!")
			self._valid_args = False
			return

		if ArgumentHandler.output_directory is not None and ArgumentHandler.project_directory is None:
			directory_out = ArgumentHandler.output_directory.rstrip('/')
			if not os.path.isdir(directory_out):
				directory_out = os.path.dirname(directory_out)

			if not os.path.isdir(directory_out):
				self._logger.error("'-od' Directory does not exist: '{}'".format(directory_out))
				self._valid_args = False
				return
			if ArgumentHandler.project_directory is None:
				if ArgumentHandler.stage < 2:
					ts = time.time()
					folder_name = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M_%S')
					ArgumentHandler.project_directory = os.path.join(ArgumentHandler.output_directory, folder_name)
					os.mkdir(ArgumentHandler.project_directory)
					self._logger.info("Project directory created: '{}'".format(ArgumentHandler.project_directory))
				else:
					# search previous project
					list_of_items = os.listdir(ArgumentHandler.output_directory)
					list_of_folder_directories = []
					for item in list_of_items:
						folder_item = os.path.join(ArgumentHandler.output_directory, item)
						if os.path.isdir(folder_item):
							mg_fasta = os.path.join(folder_item, ArgumentHandler.file_mg_16s)
							if os.path.isfile(mg_fasta):
								list_of_folder_directories.append(folder_item)

					if len(list_of_folder_directories) == 0:
						self._logger.error("'-o' No valid Project directory found at '{}'".format(directory_out))
						self._valid_args = False
						return
					elif len(list_of_folder_directories) == 1:
						ArgumentHandler.project_directory = list_of_folder_directories[0]
						self._logger.info("Selected project directory: '{}'".format(ArgumentHandler.project_directory))
					else:
						self._logger.error("'-o' Several valid project directories found at '{}'. Please specify!".format(directory_out))
						self._valid_args = False
						return


		if not os.path.isdir(ArgumentHandler.project_directory) and ArgumentHandler.stage == 0:
			os.mkdir(ArgumentHandler.project_directory)

		if ArgumentHandler.output_directory is None:
			ArgumentHandler.output_directory = os.path.dirname(ArgumentHandler.project_directory.rstrip('/'))

		if not os.path.isabs(ArgumentHandler.project_directory):
			ArgumentHandler.project_directory = os.path.expanduser(ArgumentHandler.project_directory)
			ArgumentHandler.project_directory = os.path.realpath(ArgumentHandler.project_directory)

		if not os.path.isabs(ArgumentHandler.output_directory):
			ArgumentHandler.output_directory = os.path.expanduser(ArgumentHandler.output_directory)
			ArgumentHandler.output_directory = os.path.realpath(ArgumentHandler.output_directory)

		if not os.path.isdir(ArgumentHandler.project_directory):
			self._logger.error("'-o' Directory does not exist: '{}'".format(ArgumentHandler.project_directory))
			self._valid_args = False
			return

		ArgumentHandler.silva_reference_directory = os.path.expanduser(ArgumentHandler.silva_reference_directory)
		if ArgumentHandler.silva_reference_directory is None:
			self._logger.error("'-silva' SILVA reference directory is required!")
			self._valid_args = False
			return
		if not os.path.isabs(ArgumentHandler.silva_reference_directory):
			ArgumentHandler.silva_reference_directory = os.path.realpath(ArgumentHandler.silva_reference_directory)
		if not os.path.isdir(ArgumentHandler.silva_reference_directory):
			self._logger.error("'-silva' Directory does not exist: '{}'".format(ArgumentHandler.silva_reference_directory))
			self._valid_args = False
			return

		for sub_directory in ArgumentHandler._silva_ref_files:
			file_path = os.path.join(ArgumentHandler.silva_reference_directory, sub_directory)
			if not os.path.isfile(file_path):
				self._logger.error("'-silva' File does not exist: '{}'".format(file_path))
				self._valid_args = False
				return

		if ArgumentHandler.ncbi_reference_directory is None:
			self._logger.error("'-ncbi' NCBI reference directory is required!")
			self._valid_args = False
			return

		ArgumentHandler.ncbi_reference_directory = os.path.expanduser(ArgumentHandler.ncbi_reference_directory)
		if not os.path.isabs(ArgumentHandler.ncbi_reference_directory):
			ArgumentHandler.ncbi_reference_directory = os.path.join(ArgumentHandler.pipeline_directory, ArgumentHandler.ncbi_reference_directory)

		if not os.path.isdir(ArgumentHandler.ncbi_reference_directory):
			self._logger.error("'-ncbi' Directory does not exist: '{}'".format(ArgumentHandler.ncbi_reference_directory))
			self._valid_args = False
			return

		for sub_directory in ArgumentHandler._ncbi_ref_files:
			file_path = os.path.join(ArgumentHandler.ncbi_reference_directory, sub_directory)
			if not os.path.isfile(file_path):
				self._logger.error("'-ncbi' File does not exist: '{}'".format(file_path))
				self._valid_args = False
				return

		if ArgumentHandler.stage < 2:
			if ArgumentHandler.input_reference_file is None and ArgumentHandler.input_reference_fna_file is None:
				self._logger.error("'-ir' or '-irf' Reference genome maping file is required!")
				self._valid_args = False
				return
			elif ArgumentHandler.input_reference_file is None and ArgumentHandler.input_reference_fna_file is None:
				self._logger.error("'-ir' or '-irf' Reference genome maping file is required!")
				self._valid_args = False
				return
			else:  # if not os.path.isfile(ArgumentHandler.input_reference_file) and not os.path.isfile(ArgumentHandler.input_reference_fna_file):
				file_path = ArgumentHandler.input_reference_file or ArgumentHandler.input_reference_fna_file
				if not os.path.isfile(file_path):
					self._logger.error("'-ir','-irf' File does not exist: '{}'".format(file_path))
					self._valid_args = False
					return

		if ArgumentHandler.stage > 2:
			cluster_file = os.path.join(ArgumentHandler.project_directory, ArgumentHandler.file_cluster_mg_16s)
			if not os.path.isfile(cluster_file):
				self._logger.error("Mothur file with list of clusters not found at: '{}'".format(cluster_file))
				self._valid_args = False
				return

		if ArgumentHandler.input_genomes_file is None:
			self._logger.error("'-i' Unidentified genome mapping file is required!")
			self._valid_args = False
			return
		elif not os.path.isfile(ArgumentHandler.input_genomes_file):
			self._logger.error("'-i' File does not exist: '{}'".format(ArgumentHandler.input_genomes_file))
			self._valid_args = False
			return

		if ArgumentHandler.metadata_table_in is None:
			self._logger.error("'-im' Metadata file is required!")
			self._valid_args = False
			return
		elif not os.path.isfile(ArgumentHandler.metadata_table_in):
			self._logger.error("'-im' File does not exist: '{}'".format(ArgumentHandler.metadata_table_in))
			self._valid_args = False
			return

		if ArgumentHandler.metadata_table_out is None:
			basename = os.path.basename(ArgumentHandler.metadata_table_in)
			ArgumentHandler.metadata_table_out = os.path.join(ArgumentHandler.project_directory, basename+".out.csv")
			self._logger.warning("'-om' Metadata output: '{}'".format(ArgumentHandler.metadata_table_out))
			#self._logger.error("'-om' Metadata file is required!")
			#self._valid_args = False
			#return

		if ArgumentHandler.processors is None:
			self._logger.error("'-p' A number of processors is required!")
			self._valid_args = False
			return
		elif not ArgumentHandler.processors > 0:
			self._logger.error("'-p' The number of processors must be a positive number: '{}'".format(ArgumentHandler.processors))
			self._valid_args = False
			return

		if ArgumentHandler.distance_cutoff is None:
			self._logger.error("'-th' A distance cutoff is required!")
			self._valid_args = False
			return
		elif not ArgumentHandler.distance_cutoff > 0 or not ArgumentHandler.distance_cutoff <= 1:
			self._logger.error("'-th' The distance cutoff must be between 0 and 1".format(ArgumentHandler.distance_cutoff))
			self._valid_args = False
			return

		if ArgumentHandler.otu_distance is None:
			self._logger.error("'-otu' A threshold is required!")
			self._valid_args = False
			return
		elif not ArgumentHandler.otu_distance > 0 or not ArgumentHandler.otu_distance < 1:
			self._logger.error("'-otu' The number of processors must be a positive number: '{}'".format(ArgumentHandler.otu_distance))
			self._valid_args = False
			return

		if ArgumentHandler.classification_distance_minimum is None:
			self._logger.error("'-cth' A threshold is required!")
			self._valid_args = False
			return
		elif not ArgumentHandler.classification_distance_minimum > 0 or not ArgumentHandler.classification_distance_minimum < 1:
			self._logger.error("'-cth' The number of processors must be a positive number: '{}'".format(ArgumentHandler.classification_distance_minimum))
			self._valid_args = False
			return

		expected_output_size = self._expected_output_size_in_giga_byte()
		expected_tmp_size = expected_output_size
		#if ArgumentHandler.multiple_samples:
		#	expected_tmp_size /= ArgumentHandler.number_of_samples
		directory_tmp = "/tmp"
		directory_out = ArgumentHandler.project_directory
		if not os.path.isdir(directory_out):
			directory_out = os.path.dirname(directory_out.rstrip("/"))
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
			user_input = raw_input(message+" Are you sure you want to continue? [y/n]\n>").lower()
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
		sub_folders = [ArgumentHandler.folder_name_source, ArgumentHandler.folder_name_tools]
		if ArgumentHandler.pipeline_directory is None:
			self._logger.error("Pipeline directory is required!")
			self._valid_args = False
			return
		elif not os.path.isdir(ArgumentHandler.pipeline_directory):
			self._logger.error("Pipeline directory does not exist: '{}'".format(ArgumentHandler.pipeline_directory))
			self._valid_args = False
			return

		for sub_folder in sub_folders:
			sub_dir = os.path.join(ArgumentHandler.pipeline_directory, sub_folder)
			if not os.path.isdir(sub_dir):
				self._logger.error("Bad pipeline directory: '{}'".format(ArgumentHandler.pipeline_directory))
				self._valid_args = False
				return

		#if ArgumentHandler.config_file_path is None:
		#	ArgumentHandler.config_file_path = os.path.join(ArgumentHandler.pipeline_directory, ArgumentHandler._filename_config_default)
		if not os.path.isfile(ArgumentHandler.config_file_path):
			self._logger.error("'-c' File does not exist: '{}'".format(ArgumentHandler.config_file_path))
			self._valid_args = False
			return

		self._config = Config(ArgumentHandler.config_file_path, self._logger)
		sections = ["Main", "MarkerGeneExtraction", "MarkerGeneClustering", "MarkerGeneClassification"]
		missing_section = self._config.has_missing_section(sections)
		if missing_section:
			self._logger.error("Missing section '{}' in the configuration file.".format(missing_section))
			self._valid_args = False
			return
		#if ArgumentHandler.pipeline_directory is None:
		#	ArgumentHandler.pipeline_directory = self._get_config_value("main", "pipeline_dir")

		if ArgumentHandler.output_directory is None:
			ArgumentHandler.output_directory = self._config.get_value("Main", "output_directory", verbose=False)

		if ArgumentHandler.project_directory is None:
			ArgumentHandler.project_directory = self._config.get_value("Main", "project_directory", verbose=False)

		if ArgumentHandler.processors is None:
			ArgumentHandler.processors = self._config.get_value("Main", "processors", is_digit=True)

		if ArgumentHandler.input_reference_file is None:
			ArgumentHandler.input_reference_file = self._config.get_value("MarkerGeneExtraction", "input_reference_file")

		if ArgumentHandler.input_reference_file is None:
			ArgumentHandler.input_reference_file = self._config.get_value("MarkerGeneExtraction", "input_reference_file")

		if ArgumentHandler.hmmer is None:
			ArgumentHandler.hmmer = self._config.get_value("MarkerGeneExtraction", "hmmer", is_digit=True)

		if ArgumentHandler.input_genomes_file is None:
			ArgumentHandler.input_genomes_file = self._config.get_value("MarkerGeneExtraction", "input_genomes_file")

		if ArgumentHandler.binary_mothur is None:
			ArgumentHandler.binary_mothur = self._config.get_value("MarkerGeneClustering", "mothur")

		if ArgumentHandler.metadata_table_in is None:
			ArgumentHandler.metadata_table_in = self._config.get_value("MarkerGeneClustering", "metadata_table_in")

		if ArgumentHandler.metadata_table_out is None:
			ArgumentHandler.metadata_table_out = self._config.get_value("MarkerGeneClustering", "metadata_table_out")

		if ArgumentHandler.silva_reference_directory is None:
			ArgumentHandler.silva_reference_directory = self._config.get_value("MarkerGeneClustering", "silva_reference_directory")
		else:
			self._logger.error("Not getting SILVA config for stupid reason!!")

		if ArgumentHandler.ncbi_reference_directory is None:
			ArgumentHandler.ncbi_reference_directory = self._config.get_value("MarkerGeneClassification", "ncbi_reference_directory")

		#if ArgumentHandler.genome_sample_size is None or not ArgumentHandler.genome_sample_size > 0:
		#	ArgumentHandler.genome_sample_size = self._get_config_value("sample", "num_genomes", True)

	@staticmethod
	def _free_space_in_giga_bytes(directory="/tmp"):
		if not os.path.isdir(directory):
			return 0
		statvfs = os.statvfs(directory)
		free_space = statvfs.f_frsize * statvfs.f_bfree
		return free_space / float(1024*1024*1024)

	@staticmethod
	def _expected_output_size_in_giga_byte():
		expected_output_size = 0
		return expected_output_size

	def _read_options(self, options):
		config_file = options.config_file
		if config_file is not None:
			if not os.path.isabs(config_file):
				#config_file = os.path.join(ArgumentHandler.pipeline_directory, ArgumentHandler.folder_name_tools, config_file)
				config_file = os.path.realpath(config_file)
			if not os.path.isfile(config_file):
				self._logger.error("File does not exist: '{}'".format(options.config_file))
				self._valid_args = False
				return
		ArgumentHandler.config_file_path = config_file
		ArgumentHandler._verbose = options.verbose
		ArgumentHandler._debug_mode = options.debug_mode
		ArgumentHandler._logging = options.logging
		ArgumentHandler.stage = options.stage
		ArgumentHandler.processors = options.processors

		ArgumentHandler.input_reference_file = options.input_reference_file
		ArgumentHandler.input_reference_fna_file = options.input_reference_fna_file
		ArgumentHandler.input_genomes_file = options.input_genomes
		ArgumentHandler.output_directory = options.output_directory
		ArgumentHandler.project_directory = options.project_directory
		ArgumentHandler.metadata_table_in = options.metadata_table_in
		ArgumentHandler.metadata_table_out = options.metadata_table_out
		ArgumentHandler.distance_cutoff = options.threshold
		ArgumentHandler.otu_distance = options.otu_distance
		ArgumentHandler.classification_distance_minimum = options.classification_distance
		#ArgumentHandler.ncbi_reference_directory = options.ncbi_reference_directory
		#ArgumentHandler.silva_reference_directory = options.silva_reference_directory
		#ArgumentHandler. = options.

	def _get_parser_options(self, args=None):
		description = "Pipeline for the extraction of marker genes, clustering and taxonomic classification"
		self.parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
		self.parser.add_argument("-verbose", "--verbose", action='store_true', default=False, help="display more information!")
		self.parser.add_argument("-debug", "--debug_mode", action='store_true', default=False, help="activate DEBUG modus. tmp folders will not be deleted!")
		self.parser.add_argument("-log", "--logging", action='store_true', default=False, help="pipeline output will written to a log file")
		#self.parser.add_argument("-qc", "--quality_check", action='store_true', default=False, help="")
		self.parser.add_argument("-p", "--processors", default=None, type=int,
							help="number of processors to be used. >40 recommended.")
		self.parser.add_argument("-s", "--stage", default=0, type=int, choices=[0, 1, 2, 3, 4], help='''available options: 0-4:
0 -> Full run through,
1 -> Marker gene extraction,
2 -> Gene alignment and clustering,
3 -> Classification of Genomes and novelty prediction
4 -> Average Nucleotide Identity calculation
Default: 0
''')
		#group = parser.add_mutually_exclusive_group()
		group_input = self.parser.add_argument_group("input/output")
		group_input.add_argument("-ir", "--input_reference_file", default=None, type=str,
							help="""path to a file containing list of reference genomes
Format: <genome_id>\\t<path>
No column names!""")
		group_input.add_argument("-irf", "--input_reference_fna_file", default=None, type=str,
							help="path to a fasta file containing the 16S marker genes of the reference genomes")
		group_input.add_argument("-i", "--input_genomes", default=None, type=str,
							help="""path to a file containing list of unidentified genomes
Format: <genome_id>\\t<path>
No column names!""")
		group_input.add_argument("-c", "--config_file", type=str, default=None, help="path to the configuration file of the pipeline")
		group_out = group_input.add_mutually_exclusive_group()
		group_out.add_argument("-od", "--output_directory", default=None, type=str,
							help="folder in which a subfolder will created for the output")
		group_out.add_argument("-o", "--project_directory", default=None, type=str,
							help="directory containing found marker genes and also a file in mothur format containing the clustering")
		#group_input.add_argument("-sivla", "--silva_reference_directory", default=None, type=str,
		#					help="Directory that contains the SILVA reference files, alignment, distance-matrix and name file")
		#group_input.add_argument("-ncbi", "--ncbi_reference_directory", default=None, type=str,
		#					help="Directory that contains the NCBI taxonomy dump")
		group_input.add_argument("-im", "--metadata_table_in", default=None, type=str,
							help="path to file containing tab separated list of unidentified genomes")
		group_input.add_argument("-om", "--metadata_table_out", default=None, type=str,
							help="path to file containing tab separated list of genomes and their file path")

		group_clustering = self.parser.add_argument_group("clustering")
		group_clustering.add_argument("-th", "--threshold", default=0.04, type=float,
							help="only distances up to the threshold will be calculated. Default: 0.04")
		group_clustering.add_argument("-otu", "--otu_distance", default=0.03, type=float,
							help="genetic distances at which cluster will be used as otus. Default: 0.03")
		group_clustering.add_argument("-cth", "--classification_distance", default=0.02, type=float,
							help="minimum distance for classification. Default: 0.02")

		if args is None:
			return self.parser.parse_args()
		else:
			return self.parser.parse_args(args)
