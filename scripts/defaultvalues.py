__author__ = 'hofmann'

import os
import tempfile
import random
import scripts
from scripts.Validator.validator import Validator, DefaultLogging
from scripts.ComunityDesign.communitydesign import Community


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
        pipeline_dir = os.path.dirname(self._validator.get_full_path(os.path.dirname(scripts.__file__)))
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
        self._DEFAULT_sample_size_in_base_pairs = 1 * 1000000000

        self._DEFAULT_read_simulator_type = 'art'
        self._DEFAULT_error_profile = 'hi150'
        self._DEFAULT_fragment_size_standard_deviation_in_bp = 27
        self._DEFAULT_fragments_size_mean_in_bp = 270

        # ############
        # [sampledesign]
        # ############
        self._DEFAULT_strain_simulation_template = os.path.join(
            pipeline_dir, 'scripts', 'StrainSimulationWrapper', 'sgEvolver', 'simulation_dir')
        self._DEFAULT_number_of_samples = 2
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
