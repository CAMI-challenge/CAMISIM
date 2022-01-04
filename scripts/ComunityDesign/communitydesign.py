__author__ = 'hofmann'
__version__ = '0.0.6'

import os
import random
import numpy.random as np_random
import tempfile
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.StrainSimulationWrapper.strainsimulationwrapper import StrainSimulationWrapper
from scripts.StrainSelector.strainselector import StrainSelector
from scripts.PopulationDistribution.populationdistribution import PopulationDistribution
from scripts.GenomePreparation.genomepreparation import GenomePreparation
from scripts.Validator.validator import Validator


# ##################################
#
#          Community
#
# ##################################


class Community(Validator):

    def __init__(
            self, identifier, genomes_total, genomes_real, limit_per_otu, file_path_metadata_table,
            file_path_genome_locations, file_path_gff_locations, ratio, mode,
            log_mu, log_sigma, gauss_mu=None, gauss_sigma=None,
            logfile=None, verbose=True, debug=False):
        """
        Accumulation of all community related information

        @param identifier: Community identifier
        @type identifier: str | unicode
        @param genomes_total: Total amount of genomes to be drawn from this community
        @type genomes_total: int
        @param genomes_real: Amount of real genomes to be drawn, rest will drawn from simulated ones
        @type genomes_real: int
        @param limit_per_otu: A Maximum for drawn genomes belonging to the same otu, unless more are required to be drawn
        @type limit_per_otu: int
        @param file_path_metadata_table: Table of Metadata for each genome of the community
        @type file_path_metadata_table: str | unicode
        @param file_path_genome_locations: Format: 'id \t file path to fasta file'
        @type file_path_genome_locations: str | unicode
        @param file_path_gff_locations: Format: 'id \t file path to gff file'
        @type file_path_gff_locations: str | unicode
        @param ratio: If one comm. has ratio=1 and another has ration=2, the other community will be twice the size
        @type ratio: int | long | float
        @param mode: Valid: 'replicates', 'timeseries_normal', 'timeseries_lognormal', 'differential'
        @type mode: str | unicode
        @param log_mu: Mean of drawn log distribution
        @type log_mu: int | long | float
        @param log_sigma: Standard deviation of log distribution
        @type log_sigma: int | long | float
        @param gauss_mu: Mean of drawn gauss distribution
        @type gauss_mu: int | long | float
        @param gauss_sigma: Standard deviation of gauss distribution
        @type gauss_sigma: int | long | float
        @param logfile: file handler or file path to a log file
        @type logfile: file | FileIO | StringIO | str
        @param verbose: More output and user interaction is enabled.
        @type verbose: bool
        @param debug: Display debug messages
        @type debug: bool
        """
        assert genomes_real is None or genomes_real <= genomes_total
        assert mode is None or mode in PopulationDistribution.get_valid_modes()
        if verbose is None:
            verbose = False
        super(Community, self).__init__(label="Community", logfile=logfile, verbose=verbose, debug=debug)

        if genomes_real is None:
            genomes_real = genomes_total
        self.genomes_real = genomes_real
        self.genomes_total = genomes_total
        self.limit_per_otu = limit_per_otu
        self.file_path_metadata_table = self.get_full_path(file_path_metadata_table)
        self.file_path_genome_locations = self.get_full_path(file_path_genome_locations)
        self.file_path_gff_locations = None
        if file_path_gff_locations is not None:
            self.file_path_gff_locations = self.get_full_path(file_path_gff_locations)
        self.ratio = ratio
        self.log_mu = log_mu
        self.log_sigma = log_sigma
        self.gauss_mu = gauss_mu
        self.gauss_sigma = gauss_sigma
        self.mode = mode
        self.simulate_strains = False
        if genomes_real and genomes_real < genomes_total:
            self.simulate_strains = True
        self.verbose = verbose
        self.id = identifier

    def has_valid_values(self):
        if not self.validate_characters(self.id) or self.id == '':
            return False

        if not self.validate_characters(self.mode) or self.mode == '':
            return False

        if not self.validate_number(self.genomes_total, self.genomes_real):
            return False

        if not self.validate_number(self.genomes_real, 1, self.genomes_total):
            return False

        if not self.validate_number(self.ratio, 0, zero=False):
            return False

        if not self.validate_number(self.log_mu, 0, zero=False):
            return False

        if not self.validate_number(self.log_sigma, 0, zero=True):
            return False

        if not self.validate_number(self.gauss_mu):
            return False

        if not self.validate_number(self.gauss_sigma):
            return False

        if not self.validate_number(self.limit_per_otu, 1):
            return False

        if not self.validate_file(self.file_path_metadata_table):
            return False

        if not self.validate_file(self.file_path_genome_locations):
            return False

        if self.file_path_gff_locations and not self.validate_file(self.file_path_gff_locations):
            return False

        return True

# ##################################
#
#          CommunityDesign
#
# ##################################


class CommunityDesign(GenomePreparation):
    """
        For the design of an artificial community
    """
    # _filename_distribution_comunity = "distribution_{comunity_index}_{sample_index}.txt"
    _filename_distribution_comunity_joint = "distribution_{sample_index}.txt"

    # TODO: plasmids within genome files
    # used_genomes_with_plasmids[genome_id] = random.randint(7, 10)
    # distribution = str(int(distribution) * factor)
    def __init__(
        self, column_name_genome_id="genome_ID", column_name_otu="OTU",
        column_name_novelty_category="novelty_category", column_name_ncbi="NCBI_ID", column_name_source="source",
        max_processors=1, tmp_dir=None, logfile=None, verbose=True, debug=False, seed=None):
        """
        @param column_name_genome_id: Column name of genome ids in the metadata table
        @type column_name_genome_id: str | unicode
        @param column_name_otu: Column name of otu ids in the metadata table
        @type column_name_otu: str | unicode
        @param column_name_novelty_category: Column name of novelty category in the metadata table
        @type column_name_novelty_category: str | unicode
        @param column_name_ncbi: Column name of taxonomic id assignment in the metadata table
        @type column_name_ncbi: str | unicode
        @param column_name_source: Column name of 'source' in the metadata table
        @type column_name_source: str | unicode
        @param max_processors: maximum number of processors available to be used
        @type max_processors: long | int
        @param tmp_dir: working directory or place temporary files can be stored
        @type tmp_dir: str | unicode
        @param logfile: file handler or file path to a log file
        @type logfile: file | FileIO | StringIO | str
        @param verbose: Not verbose means that only warnings and errors will be past to stream
        @type verbose: bool
        @param debug: Display debug messages
        @type debug: bool
        """
        super(CommunityDesign, self).__init__(label="CommunityDesign", logfile=logfile, verbose=verbose, debug=debug)
        if seed is not None:
            random.seed(seed)
            np_random.seed(abs(hash(seed)) % 4294967295)  # numpy accepts only 32 bit integers

        # self._seed = seed
        # self._filename_distribution = filename_prefix_distribution + "{index}.txt"
        self._column_name_genome_id = column_name_genome_id
        self._column_name_otu = column_name_otu
        self._column_name_novelty_category = column_name_novelty_category
        self._column_name_source = column_name_source
        self._column_name_ncbi = column_name_ncbi

        assert isinstance(max_processors, int)
        assert max_processors > 0
        self._max_processors = max_processors

        if tmp_dir is None:
            tmp_dir = tempfile.gettempdir()
        self._tmp_dir = tmp_dir
        assert self.validate_dir(self._tmp_dir)

    @staticmethod
    def get_distribution_file_paths(directory, number_of_samples):
        """
        Generate directory paths for each sample

        @param directory: Output stream
        @type directory: str | unicode
        @param number_of_samples: Number of samples
        @type number_of_samples: int | long

        @return: list of directories
        @rtype: list[str | unicode]
        """
        file_path = os.path.join(directory, CommunityDesign._filename_distribution_comunity_joint)
        return [file_path.format(sample_index=sample_index) for sample_index in range(number_of_samples)]

    @staticmethod
    def _write_distribution_file(stream_out, genome_id_to_abundance):
        """
        Write abundance file for each sample

        @param stream_out: Output stream
        @type stream_out: file | FileIO | StringIO
        @param genome_id_to_abundance: Drawn distribution for each genome id
        @type genome_id_to_abundance: dict[str|unicode, list[float]]
        """
        for genome_id in genome_id_to_abundance:
            distributions = [str(abundance) for abundance in genome_id_to_abundance[genome_id]]
            stream_out.write("{id}\t{distr}\n".format(id=genome_id, distr='\t'.join(distributions)))

    def design_community(
        self, file_path_distributions, community, number_of_samples, metadata_table,
        directory_out_metadata, directory_in_template=None):
        """
        Design artificial community, of a specific design, with different distributions for each sample

        @param file_path_distributions: File path where distributions will be written to
        @type file_path_distributions: str | unicode
        @param community: Input data for the creation of a community
        @type community: Community
        @param number_of_samples: Amount of samples to be simulated
        @type number_of_samples: int
        @param metadata_table: Will contain metadata of all (simulated) genomes/plasmids drawn
        @type metadata_table: MetadataTable
        @param directory_out_metadata: Metadata tables of separated by chosen and not chosen genomes are written to here
        @type directory_out_metadata: str | unicode
        @param directory_in_template: contains template data for strain simulation
        @type directory_in_template: str | unicode

        @return: Dictionary with drawn genome ids as key and file paths as value
        @rtype: dict[str|unicode, str|unicode]
        """
        assert isinstance(community, Community)
        assert isinstance(metadata_table, MetadataTable)

        number_of_strains = community.genomes_total

        # pick how much a strain will be simulated
        genome_amounts = []
        strain_simulation = None
        if community.simulate_strains:
            strain_simulation = StrainSimulationWrapper(
                executable_sim=None,
                directory_template=directory_in_template,
                column_name_gid=self._column_name_genome_id,
                column_name_ncbi=self._column_name_ncbi,
                column_name_source=self._column_name_source,
                separator='\t',
                filename_prefix="simulated_",
                keep_original=True,
                max_processors=self._max_processors,
                tmp_dir=self._tmp_dir,
                logfile=self._logfile, verbose=self._verbose, debug=self._debug,
                # seed=self._seed
                )

            probability = None  # 1-options.communities[community_id]["evolve"]
            genome_amounts = strain_simulation.get_genome_amounts(
                probability=probability,
                max_genome_amount=community.genomes_total,
                num_real_genomes=community.genomes_real,
                silent=not community.verbose
            )
            number_of_strains = len(genome_amounts)

        # draw strains
        self._logger.info("Drawing strains.")
        metadata_table_community = MetadataTable(logfile=self._logfile, verbose=self._verbose)
        metadata_table_community.read(community.file_path_metadata_table, column_names=True)
        strain_selector = StrainSelector(
            column_name_genome_id=self._column_name_genome_id,
            column_name_otu=self._column_name_otu,
            column_name_novelty_category=self._column_name_novelty_category,
            logfile=self._logfile, verbose=self._verbose, debug=self._debug
            )
        list_of_drawn_genome_id = strain_selector.get_drawn_genome_id(
            metadata_table=metadata_table_community,
            number_of_strains=number_of_strains,
            number_of_strains_per_otu=community.limit_per_otu
            )

        # write unused data to separate file
        old_base_name = os.path.basename(community.file_path_metadata_table)
        file_prefix, extention = os.path.splitext(old_base_name)
        new_file_name = "unused_c{index}_{prefix}{ext}".format(
            prefix=file_prefix,
            index=community.id,
            ext=extention)
        metadata_new_file_path = os.path.join(directory_out_metadata, new_file_name)
        metadata_table_community.write(
            metadata_new_file_path,
            exclude=True,
            value_list=list_of_drawn_genome_id,
            key_column_name=self._column_name_genome_id,
            column_names=True)

        # get path for every genome
        genome_id_to_file_path_gff = None
        if community.file_path_gff_locations:
            genome_id_to_file_path_gff = self._get_genome_id_to_path_map(
                community.file_path_gff_locations, list_of_drawn_genome_id)
        genome_id_to_path_map = self._get_genome_id_to_path_map(
            community.file_path_genome_locations, list_of_drawn_genome_id)

        # concatenate
        metadata_table_community.reduce_rows_to_subset(list_of_drawn_genome_id, self._column_name_genome_id)
        metadata_table.concatenate(metadata_table_community, strict=False)

        # validate correct format of files
        self._logger.info("Validating raw sequence files!")
        assert self.validate_format(
            list_of_file_paths=genome_id_to_path_map.values(),
            file_format="fasta",
            sequence_type="dna",
            ambiguous=True
            ), "Validation of file format failed!"

        # simulate diversity around strains
        if community.simulate_strains:
            genome_id_to_amounts = strain_simulation.get_genome_id_to_amounts(list_of_drawn_genome_id, genome_amounts)
            strain_simulation.simulate_strains(
                meta_table=metadata_table,
                genome_id_to_amounts=genome_id_to_amounts,
                genome_id_to_file_path_genome=genome_id_to_path_map,
                genome_id_to_file_path_gff=genome_id_to_file_path_gff)
            # adopt new list that includes simulated strains
            self._logger.info("Validating simulated sequence files!")
            for genome_id, file_path in genome_id_to_path_map.items():
                if genome_id in list_of_drawn_genome_id:
                    continue
                assert self.validate_sequence_file(
                    file_path,
                    file_format="fasta",
                    sequence_type="dna",
                    ambiguous=True)
            list_of_drawn_genome_id = genome_id_to_path_map.keys()

        # get community distributions
        population_distribution = PopulationDistribution(
            logfile=self._logfile, verbose=self._verbose, debug=self._debug)
        list_of_distributions = population_distribution.get_lists_of_distributions(
            size_of_population=len(list_of_drawn_genome_id),
            number_of_samples=number_of_samples,
            modus=community.mode,
            log_mu=community.log_mu, log_sigma=community.log_sigma,
            gauss_mu=community.gauss_mu, gauss_sigma=community.gauss_sigma,
            view_distribution=community.verbose
        )

        # move and clean up files (removes sequence description)
        # genome_id_to_total_length = self.move_genome_files(
        #     genome_id_to_path_map,
        #     directory_output=directory_out_genomes,
        #     sequence_min_length=min_sequence_length,
        #     set_of_sequence_names=set_of_sequence_names)

        # write distribution file
        # genome_id_to_distributions = self._get_genome_id_to_distributions(list_of_drawn_genome_id, list_of_distributions)
        assert len(list_of_drawn_genome_id) == len(list_of_distributions)
        genome_id_to_distributions = dict(zip(list_of_drawn_genome_id, list_of_distributions))

        # genome_id_to_file_name = self._get_genome_id_to_file_name(genome_id_to_path_map)
        with open(file_path_distributions, 'w') as stream_out:
            self._write_distribution_file(stream_out=stream_out, genome_id_to_abundance=genome_id_to_distributions)
        return genome_id_to_path_map

    @staticmethod
    def _get_genome_id_to_file_name(genome_id_to_path_map):
        """
        Extract file names for reuse in new folder

        @attention: Preserving original name might will cause issue if not unique

        @param genome_id_to_path_map: Dictionary with file path by genome ids
        @type genome_id_to_path_map: dict[str|unicode, str|unicode]

        @return: Map of genome id to a filename
        @rtype : dict[str|unicode, str|unicode]
        """
        set_of_file_names = set()
        genome_id_to_file_name = {}
        for genome_id, file_path in genome_id_to_path_map.items():
            filename = os.path.basename(file_path)
            assert filename not in set_of_file_names, "Filename '{}' is not unique!".format(set_of_file_names)
            set_of_file_names.add(filename)
            genome_id_to_file_name[genome_id] = filename
        return genome_id_to_file_name

    #####################
    #
    # merge communities
    #
    #####################

    def merge_communities(
        self, list_of_communities, list_of_comunity_distribution_file_paths, index_sample, file_path_output):
        """
        Combine distributions of communities and adjust them according to their ratio.

        @param list_of_communities: List of community inputs
        @type list_of_communities: list[Community]
        @param list_of_comunity_distribution_file_paths: List of distributions
        @type list_of_comunity_distribution_file_paths: list[str | unicode]
        @param index_sample: Index of sample
        @type index_sample: int | long
        @param file_path_output: Sample distribution file path
        @type file_path_output: str | unicode

        @return: Nothing
        @rtype: None
        """
        assert isinstance(list_of_communities, list)
        for community in list_of_communities:
            assert isinstance(community, Community)
        # assert isinstance(metadata_table, MetadataTable)

        # read communities and adapt to ratio
        list_of_community_total_abundance = [0] * len(list_of_communities)
        sample_total_abundance = 0

        genomes = set()
        metadata_table_community = MetadataTable(logfile=self._logfile, verbose=self._verbose)
        for index_community, file_path in enumerate(list_of_comunity_distribution_file_paths):
            community_distribution = metadata_table_community.parse_file(file_path, column_names=False)
            for row in community_distribution:
                genome_id = row[0]
                if genome_id in genomes:
                    raise ValueError("Genome id '{}' not unique".format(genome_id))
                genomes.add(genome_id)
                abundance = row[index_sample+1]
                list_of_community_total_abundance[index_community] += float(abundance)  # * float(sequence_info[4])
            community_distribution.close()

        for index_community, _ in enumerate(list_of_comunity_distribution_file_paths):
            sample_total_abundance += list_of_community_total_abundance[index_community]

        # out.append(read_communities[0][0])
        list_of_community_factor = [0.0] * len(list_of_communities)

        for index_community, _ in enumerate(list_of_comunity_distribution_file_paths):
            ratio = float(list_of_communities[index_community].ratio)
            community_total_abundance = float(list_of_community_total_abundance[index_community])
            current_proportion_in_sample = community_total_abundance / float(sample_total_abundance)
            list_of_community_factor[index_community] = ratio / current_proportion_in_sample
            # self.update_community(communities[index_community], factor)

        # join communities
        communities = []
        for index_community, file_path in enumerate(list_of_comunity_distribution_file_paths):
            communities.append(metadata_table_community.parse_file(file_path, column_names=False))

        # print_ratios(communities)
        with open(file_path_output, 'w') as stream_output:
            self._write_joined_community(
                communities,
                list_of_community_factor,
                index_sample,
                stream_output)

    @staticmethod
    def _write_joined_community(communities, list_of_community_factor, index_sample, stream_output):
        """
        Stream out joined distribution for a sample

        @param communities: list of iterators over row cells from community abundance files
        @type communities: list[generator]
        @param list_of_community_factor: multiplication factor for each community to get right ratio
        @type list_of_community_factor: list[float]
        @param index_sample: Index of sample
        @type index_sample: int | long
        @param stream_output: joined distribution information output
        @type stream_output: file | FileIO | StringIO | str
        """
        gid_to_abundance = {}
        total_abundance = 0.0
        line_format = "{gid}\t{abundance}\n"
        for community_index, community in enumerate(communities):
            factor = list_of_community_factor[community_index]
            for row in community:
                genome_id = row[0]
                abundance = row[index_sample+1]
                gid_to_abundance[genome_id] = float(abundance) * factor
                total_abundance += gid_to_abundance[genome_id]
            community.close()
        for genome_id, abundance in gid_to_abundance.items():
            stream_output.write(line_format.format(
                gid=genome_id,
                abundance=float(abundance) / total_abundance  # saving relative abundance
            ))

    def design_samples(
        self, list_of_communities, metadata_table, list_of_file_paths_distribution, directory_out_metadata,
        directory_in_template=None):
        """
        Design artificial community, of a specific design, with different distributions for each sample

        @param list_of_file_paths_distribution: Output file path list for distribution files
        @type list_of_file_paths_distribution: list[str | unicode]
        @param list_of_communities: List of input data for the creation of a community
        @type list_of_communities: list[Community]
        @param metadata_table: Will contain metadata of all (simulated) genomes/plasmids drawn
        @type metadata_table: MetadataTable
        @param directory_out_metadata: Metadata tables of separated by chosen and not chosen genomes are written to here
        @type directory_out_metadata: str | unicode
        @param directory_in_template: contains template data for strain simulation
        @type directory_in_template: str | unicode

        @return: Dictionary with drawn genome ids as key and file paths as value
        @rtype: tuple[dict[str|unicode, str|unicode]
        """
        assert isinstance(list_of_communities, list)
        assert isinstance(list_of_file_paths_distribution, list)
        for community in list_of_communities:
            assert isinstance(community, Community)
        assert isinstance(metadata_table, MetadataTable)

        list_of_comunity_distribution_file_paths = []

        merged_genome_id_to_path_map = {}
        for community in list_of_communities:
            file_path_output_comunity = tempfile.mktemp(dir=self._tmp_dir)  # insecure
            list_of_comunity_distribution_file_paths.append(file_path_output_comunity)
            genome_id_to_path_map = self.design_community(
                file_path_distributions=file_path_output_comunity,
                community=community,
                number_of_samples=len(list_of_file_paths_distribution),
                metadata_table=metadata_table,
                directory_out_metadata=directory_out_metadata,
                directory_in_template=directory_in_template)
            merged_genome_id_to_path_map.update(genome_id_to_path_map)

        for index_sample, file_path_output in enumerate(list_of_file_paths_distribution):
            self.merge_communities(
                list_of_communities, list_of_comunity_distribution_file_paths, index_sample, file_path_output)

        # delete now obsolete files
        if not self._debug:
            for file_path in list_of_comunity_distribution_file_paths:
                if os.path.isfile(file_path):
                    os.remove(file_path)

        return merged_genome_id_to_path_map
