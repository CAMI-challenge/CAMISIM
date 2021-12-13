__author__ = 'hofmann'

import sys
from scripts.configparserwrapper import ConfigParserWrapper
from scripts.Validator.validator import Validator
from scripts.ComunityDesign.communitydesign import Community
from scripts.defaultvalues import DefaultValues


class ConfigFileHandler(DefaultValues):
    """
    Reading and writing config file

    @type _list_of_communities: list[Community]
    """
    # internal variables not set in config
    _file_name_config = "config.ini"
    _ncbi_ref_files = ["nodes.dmp", "merged.dmp", "names.dmp"]

    def __init__(self, label="ConfigFileHandler", logfile=None, verbose=False, debug=False):
        super(ConfigFileHandler, self).__init__(label=label, logfile=logfile, verbose=verbose, debug=debug)
        self._validator = Validator(logfile=logfile, verbose=verbose, debug=debug)

    def _read_config(self, file_path_config):
        """
        Read parameter from configuration file.

        @rtype: bool
        """
        # TODO: check that all keys options make sense
        self._config = ConfigParserWrapper(logfile=self._logfile, verbose=self._verbose)
        if not self._validator.validate_file(file_path_config, key="Configuration file"):
            self._valid_args = False
            return
        self._config.read(file_path_config)

        # ##########
        # [Main]
        # ##########

        section = None  # "Main"
        if self._phase is None:
            self._phase = self._config.get_value("phase", is_digit=True, silent=True)

        if self._seed is None:
            self._seed = self._config.get_value("seed", silent=True)

        if self._max_processors is None:
            self._max_processors = self._config.get_value("max_processors", is_digit=True, silent=True)

        if self._dataset_id is None:
            self._dataset_id = self._config.get_value("dataset_id", silent=True)

        if self._directory_output is None:
            self._directory_output = self._config.get_value("output_directory", is_path=True)

        if self._tmp_dir is None:
            config_value = self._config.get_value("temp_directory", is_path=True, silent=True)
            if config_value is not None:
                assert self._validator.validate_dir(config_value)
                self._tmp_dir = config_value

        self._phase_gsa = self._config.get_value("gsa", is_boolean=True, silent=True)
        self._phase_pooled_gsa = self._config.get_value("pooled_gsa", is_boolean=True, silent=True)

        self._compresslevel = self._config.get_value("compress", is_digit=True, silent=True)

        self._phase_anonymize = self._config.get_value("anonymous", is_boolean=True, silent=True)

        # ##########
        # [ReadSimulator]
        # ##########

        section = None  # "ReadSimulator"
        if self._sample_size_in_base_pairs is None:
            config_value = self._config.get_value("size", is_digit=True, silent=True)
            if config_value is not None:
                self._sample_size_in_base_pairs = config_value * self._base_pairs_multiplication_factor

        if self._read_simulator_type is None:
            self._read_simulator_type = self._config.get_value("type", silent=True)

        if self._executable_samtools is None:
            self._executable_samtools = self._config.get_value("samtools", is_path=True, silent=True)

        if self._executable_readsim is None:
            self._executable_readsim = self._config.get_value("readsim", silent=True, is_path=True)

        if self._directory_error_profiles is None:
            self._directory_error_profiles = self._config.get_value("error_profiles", silent=True, is_path=True)

        if self._error_profile is None:
            self._error_profile = self._config.get_value("profile", silent=True)
            
        if self._custom_profile_filename is None:
            self._custom_profile_filename = self._config.get_value("base_profile_name", silent=True)
        
        if self._custom_readlength is None:
            self._custom_readlength = self._config.get_value("profile_read_length", is_digit=True, silent=True)

        if self._fragment_size_standard_deviation_in_bp is None:
            self._fragment_size_standard_deviation_in_bp = self._config.get_value(
                "fragment_size_standard_deviation", is_digit=True, silent=True)

        if self._fragments_size_mean_in_bp is None:
            self._fragments_size_mean_in_bp = self._config.get_value("fragments_size_mean", is_digit=True, silent=True)

        # ##########
        # [CommunityDesign]
        # ##########

        if self._input_list_of_file_paths_distributions is None:
            input_list_of_file_paths_distributions = self._config.get_value("distribution_file_paths", is_path=True, silent=True)
            if input_list_of_file_paths_distributions is not None:
                self._input_list_of_file_paths_distributions = input_list_of_file_paths_distributions.split(',')

        section = None  # "CommunityDesign"
        if self._directory_ncbi_taxdump is None:
            self._directory_ncbi_taxdump = self._config.get_value("ncbi_taxdump", is_path=True, silent=True)

        if self._strain_simulation_template is None:
            self._strain_simulation_template = self._config.get_value(
                "strain_simulation_template", is_path=True, silent=True)

        if self._number_of_samples is None:
            self._number_of_samples = self._config.get_value("number_of_samples", is_digit=True, silent=True)

        # if self._number_of_communities is None:
        #     self._number_of_communities = self._config.get_value('number_of_communities', is_digit=True)
        #
        # if self._number_of_communities is None:
        #     self._logger.error("Bad number of communities!")
        #     self._valid_arguments = False
        #     return

        community_sections = set()
        community_key_options = {
            "genomes_total", 'num_real_genomes', 'max_strains_per_otu', 'ratio',
            'log_mu', 'log_sigma', 'gauss_mu', 'gauss_sigma'}
        for key_options in community_key_options:
            community_sections = community_sections.union(self._config.search_sections_of(key_options))

        self._list_of_communities = []
        is_valid = True
        for community_section in community_sections:
            file_path_metadata_table = self._config.get_value('metadata', community_section, is_path=True)
            file_path_genome_locations = self._config.get_value('id_to_genome_file', community_section, is_path=True)
            file_path_gff_locations = self._config.get_value('id_to_gff_file', community_section, is_path=True, silent=True)
            mode = self._config.get_value('mode', community_section, silent=True)
            if not isinstance(file_path_metadata_table, str):
                is_valid = False
            if not isinstance(file_path_genome_locations, str):
                is_valid = False
            # if not isinstance(file_path_gff_locations, str):
            #     is_valid = False
            # if not isinstance(mode, str):
            #     is_valid = False

            if not is_valid:
                continue
            assert isinstance(file_path_metadata_table, str)
            assert isinstance(file_path_genome_locations, str)
            assert file_path_gff_locations is None or isinstance(file_path_gff_locations, str)
            assert mode is None or isinstance(mode, str)
            new_community = Community(
                identifier=community_section,
                genomes_total=self._config.get_value('genomes_total', community_section, is_digit=True),
                genomes_real=self._config.get_value('num_real_genomes', community_section, is_digit=True, silent=True),
                limit_per_otu=self._config.get_value('max_strains_per_otu', community_section, is_digit=True, silent=True),
                file_path_metadata_table=file_path_metadata_table,
                file_path_genome_locations=file_path_genome_locations,
                file_path_gff_locations=file_path_gff_locations,
                ratio=self._config.get_value('ratio', community_section, is_digit=True, silent=True),
                mode=mode,
                log_mu=self._config.get_value('log_mu', community_section, is_digit=True, silent=True),
                log_sigma=self._config.get_value('log_sigma', community_section, is_digit=True, silent=True),
                gauss_mu=self._config.get_value('gauss_mu', community_section, is_digit=True, silent=True),
                gauss_sigma=self._config.get_value('gauss_sigma', community_section, is_digit=True, silent=True),
                verbose=self._config.get_value('view', community_section, is_boolean=True, silent=True)
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
        output_stream.write("readsim={}\n".format(self._executable_readsim))
        output_stream.write("error_profiles={}\n".format(self._directory_error_profiles or ""))
        output_stream.write("samtools={}\n".format(self._executable_samtools))
        output_stream.write("profile={}\n".format(self._error_profile))
        output_stream.write("base_profile_name={}\n".format(self._custom_profile_filename or ""))
        output_stream.write("profile_read_length={}\n".format(self._custom_readlength or ""))
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
            output_stream.write("num_real_genomes={}\n".format(community.genomes_real))
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
