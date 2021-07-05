__author__ = 'hofmann'
__version__ = '0.0.3'

import os
import sys
import argparse
import tempfile
import random
import numpy.random as np_random
from scripts.projectfilefolderhandle import ProjectFileFolderHandle
from scripts.configfilehandler import ConfigFileHandler
from scripts.Archive.archive import Archive


class ArgumentHandler(ConfigFileHandler):
    """Reading pipeline configuration from file and from passed arguments"""

    _separator = None
    _file_path_config = None

    _column_name_genome_id = "genome_ID",
    _column_name_otu = "OTU",
    _column_name_novelty_category = "novelty_category",
    _column_name_ncbi = "NCBI_ID",
    _column_name_source = "source",

    def __init__(
        self, args=None, separator="\t",
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

        options = self._get_parser_options(args)
        logfile = options.logfile
        if logfile is not None:
            logfile = self._validator.get_full_path(logfile)
        super(ArgumentHandler, self).__init__(label="MetagenomeSimulationPipeline", logfile=logfile)
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
        assert isinstance(tmp_dir, str)
        assert isinstance(directory_output, str)

        if self._seed is not None:
            random.seed(self._seed)
            np_random.seed(abs(hash(self._seed)) % 4294967295)  # numpy accepts only 32 bit integers

        assert isinstance(self._directory_output, str)
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
        #     expected_output_size = 0.001
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
readsim={art_illumina}

# Directory containing error profiles for ART_Illumina
error_profiles={error_profiles}

# Supported profiles: "mi": EmpMiSeq250R, "hi": EmpHiSeq2kR, "hi150": HiSeq2500L150R
profile={error}

# Simulate samples of this size (giga base pairs)
size={gbps}

# Read simulator type (only ART is currently supported)
type={readsim_type}

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
            readsim=self._executable_readsim,
            error_profiles=self._directory_error_profiles,
            error=self._error_profile,
            gbps=float(self._sample_size_in_base_pairs)/self._base_pairs_multiplication_factor,
            readsim_type=self._read_simulator_type,
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
            self._logger.error("NCBI taxdump is required!")
            self._valid_arguments = False
        # elif not self._validator.validate_dir(
        #   self._directory_ncbi_taxdump, file_names=self._ncbi_ref_files, key='-ncbi'):
        elif not os.path.exists(self._directory_ncbi_taxdump):
            self._valid_arguments = False
        else:
            if os.path.isdir(self._directory_ncbi_taxdump):
                if not self._validator.validate_dir(self._directory_ncbi_taxdump, file_names=self._ncbi_ref_files, key='-ncbi'):
                    self._valid_arguments = False
            else:
                if not Archive().is_archive(self._directory_ncbi_taxdump):
                    self._logger.error("NCBI taxdump must be directory or compressed file! Got:\n\t'{}'".format(self._directory_ncbi_taxdump))
                    self._valid_arguments = False
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
        elif self._read_simulator_type in self._valid_read_simulators:
            if self._directory_error_profiles is None and self._read_simulator_type != 'wgsim': # wgsim does not need an error_profile folder
                self._logger.error("For %s an error profile directory is required!" % self._read_simulator_type)
                self._valid_arguments = False
            elif not self._validator.validate_dir(self._directory_error_profiles) and self._read_simulator_type != 'wgsim':
                self._valid_arguments = False
            elif self._read_simulator_type != 'wgsim':
                self._directory_error_profiles = self._validator.get_full_path(self._directory_error_profiles)

            if self._executable_readsim is None:
                self._logger.error("Read simulator executable is required!")
                self._valid_arguments = False
            elif not self._validator.validate_file(self._executable_readsim, executable=True):
                self._valid_arguments = False
            else:
                self._executable_readsim = self._validator.get_full_path(self._executable_readsim)

            if self._error_profile is None:
                self._logger.error("No error profile is given, check manual for possible choices!")
                self._valid_arguments = False

            if self._fragments_size_mean_in_bp is None and 'nanosim' not in self._read_simulator_type:
                self._logger.error("'-fmean' For the simulation with 'art', 'pbsim' and 'wgsim' the fragment size mean is required!")
                self._valid_arguments = False
            elif not self._validator.validate_number(self._fragments_size_mean_in_bp, minimum=1, key='-fmean') and 'nanosim' not in self._read_simulator_type:
                self._valid_arguments = False

            if self._fragment_size_standard_deviation_in_bp is None and self._read_simulator_type != 'nanosim':
                self._logger.error("'-fsd' For the simulation with 'art', 'pbsim' and 'wgsim' a standard_deviation of the fragments size is required!")
                self._valid_arguments = False
            elif not self._validator.validate_number(self._fragment_size_standard_deviation_in_bp, minimum=1, key='-fsd') and 'nanosim' not in self._read_simulator_type:
                self._logger.error(
                    "'-fsd' The standard_deviation of the fragments size must be a positive number: '{}'".format(
                        self._fragment_size_standard_deviation_in_bp))
                self._valid_arguments = False
        else:
            self._logger.error("The chosen read simulator %s is not supported, must be one of %s" % (self._read_simulator_type, self._valid_read_simulators))
        
        expected_output_size = self._expected_output_size_in_giga_byte()
        expected_tmp_size = expected_output_size / self._number_of_samples
        assert isinstance(self._directory_output, str)
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
            self._phase_simulate_reads = False

        if self._phase == 2:
            self._phase_validate_raw_genomes = False
            self._phase_design_community = False
            self._phase_move_and_clean_genomes = False
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
        self._verbose = not options.silent
        self._debug = options.debug_mode
        self._phase = options.phase
        self._dataset_id = options.data_set_id
        self._max_processors = options.max_processors
        self._seed = options.seed
        # self._directory_output = options.output_directory
        # self._sample_size_in_base_pairs = options.sample_size_gbp
        # if self._sample_size_in_base_pairs is not None:
        #     self._sample_size_in_base_pairs = long(options.sample_size_gbp * self._base_pairs_multiplication_factor)
        # self.read_simulator = options.read_simulator
        # self._error_profile = options.error_profile
        # self._fragment_size_standard_deviation_in_bp = options.fragment_size_standard_deviation
        # self._fragments_size_mean_in_bp = options.fragments_size_mean
        # self.plasmid_file = options.plasmid_file
        # self._number_of_samples = options.number_of_samples
        # self._phase_pooled_gsa = options.pooled_gsa

    @staticmethod
    def _get_parser_options(args=None):
        """
        Parsing of passed arguments.

        @param args: Passed arguemnts

        @return: any
        """
        parser = argparse.ArgumentParser(
            usage="python %(prog)s configuration_file_path",
            description="""
    #######################################
    #    MetagenomeSimulationPipeline     #
    #######################################

    Pipeline for the simulation of a metagenome""",
            formatter_class=argparse.RawTextHelpFormatter)

        parser.add_argument(
            "-silent", "--silent",
            action='store_true',
            default=False,  # set None if read from config
            help="Hide unimportant Progress Messages.")
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
        #     "-o", "--output_directory",
        #     type=str,
        #     default=None,
        #     help="Directory in which the simulated data will be put into")

        # ##########
        # read simulator
        # ##########

#         group_read_simulator = parser.add_argument_group('read simulator')
#         group_read_simulator.add_argument(
#             "-rs", "--read_simulator",
#             type=str,
#             default=None,
#             choices=['art', 'pirs', 'pacbio'],
#             help='''choose read simulator:
#     'pacbio',
#     'art' (Illumina),
#     'pirs' (Illumina).
# Currently only 'art' is supported.''')
#         group_read_simulator.add_argument(
#             "-ep", "--error_profile",
#             type=str,
#             default=None,
#             choices=['mi', 'hi', 'hi150'],
#             help='''Art Illumina error profiles:
#     'mi': MiSeq 250bp,
#     'hi': HiSeq 100bp,
#     'hi150': HiSeq 150bp''')
#         group_read_simulator.add_argument(
#             "-bp", "--sample_size_gbp",
#             type=float,
#             default=None,
#             help="The approximate total size the output will have in giga base pairs")
#         group_read_simulator.add_argument(
#             "-fmean", "--fragments_size_mean",
#             type=int,
#             default=None,
#             help='''Mean size of the fragments in base pairs.''')
#         group_read_simulator.add_argument(
#             "-fsd", "--fragment_size_standard_deviation",
#             type=int,
#             default=None,
#             help="Standard deviation from the mean size of fragments in base pairs.")
#
#         # ##########
#         # sample design
#         # ##########
#
#         group_community_design = parser.add_argument_group('sample design')
#         group_community_design.add_argument(
#             "-ns", "--number_of_samples",
#             type=int,
#             default=None,
#             help='''Number of samples to be simulated''')
#         group_community_design.add_argument(
#             "-pgsa", "--pooled_gsa",
#             default=None,
#             action='store_true',
#             help="Reads of all samples are pooled and assembled")
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
