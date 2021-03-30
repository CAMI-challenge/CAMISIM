__author__ = 'hofmann'
__version__ = '0.0.2.1'

import os
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy
from scripts.Validator.validator import Validator


class TaxonomicProfile(Validator):
    """
    Constructing taxonomic profiles from files with relative abundances.
    """

    _taxonomic_profile_version = "0.9.1"

    def __init__(self, taxonomy, logfile=None, verbose=True, debug=False):
        """
        @param taxonomy: taxonomy handler
        @type taxonomy: NcbiTaxonomy
        @param logfile: file handler or file path to a log file
        @type logfile: file | FileIO | StringIO | str
        @param verbose: Not verbose means that only warnings and errors will be past to stream
        @type verbose: bool
        @param debug: Display debug messages
        @type debug: bool
        """
        super(TaxonomicProfile, self).__init__(label="TaxonomicProfile", logfile=logfile, verbose=verbose, debug=debug)
        self._ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
        assert isinstance(taxonomy, NcbiTaxonomy)
        self._taxonomy = taxonomy
        self._filename_taxonomic_profile = "taxonomic_profile_{sample_index}.txt"

    def write_taxonomic_profile_from_abundance_files(
        self, metadata_table, list_of_file_paths, directory_output, sample_id=""):
        """
        Write a taxonomic profile file for each relative abundance file

        @param metadata_table: Contains metadata of all communities
        @type metadata_table: MetadataTable
        @param list_of_file_paths: List of abundance file paths
        @type list_of_file_paths: list[str | unicode]
        @param directory_output: Profiles are written in this directory
        @type directory_output: str | unicode
        @param sample_id: Identifier of a sample
        @type sample_id: str | unicode
        """
        metadata_table_tmp = MetadataTable(logfile=self._logfile, verbose=self._verbose)
        for index_abundance, file_path in enumerate(list_of_file_paths):
            community_abundance = metadata_table_tmp.parse_file(file_path, column_names=False)
            file_path_output = os.path.join(directory_output, self._filename_taxonomic_profile.format(
                sample_index=index_abundance))
            with open(file_path_output, 'w') as stream_output:
                self.write_taxonomic_profile(
                    community_abundance,
                    stream_output,
                    metadata_table,
                    sample_id)

    def write_taxonomic_profile(self, community_abundance, stream_output, metadata_table, sample_id=""):
        """
        Stream a taxonomic profile by list of relative abundances

        @param community_abundance: list of relative abundances
        @type community_abundance: generator[ dict[int|long|str|unicode, str|unicode] ]
        @param stream_output: Output of taxonomic profile
        @type stream_output: file | FileIO | StringIO
        @param metadata_table: Contains metadata of all communities
        @type metadata_table: MetadataTable
        @param sample_id: Identifier of a sample
        @type sample_id: str | unicode
        """
        assert isinstance(metadata_table, MetadataTable)
        genome_abundance = {}
        total_abundance = 0.0

        # for community in community_abundance:
        #   all_communities += community

        for genome_id, abundance in community_abundance:
            if genome_id in genome_abundance:
                raise IOError("genome id '{}' is not unique!".format(genome_id))
            genome_abundance[genome_id] = float(abundance)  # *float(total_length)
            total_abundance += genome_abundance[genome_id]

        for key, value in genome_abundance.items():
            genome_abundance[key] = value / total_abundance

        self._stream_taxonomic_profile(stream_output, genome_abundance, metadata_table, sample_id)

    def _stream_taxonomic_profile(self, stream_output, genome_id_to_percent, metadata_table, sample_id=""):
        """
        Stream a taxonomic profile by list of percentages by genome id

        @param stream_output: Output of taxonomic profile
        @type stream_output: file | FileIO | StringIO
        @param genome_id_to_percent: Percentage for each genome id
        @type genome_id_to_percent: dict[str|unicode, float]
        @param metadata_table: Contains metadata of all communities
        @type metadata_table: MetadataTable
        @param sample_id: Identifier of a sample
        @type sample_id: str | unicode
        """
        strain_id_to_genome_id = {}
        genome_id_to_strain_id = {}

        genome_id_to_taxid = metadata_table.get_map(key_column_name="genome_ID", value_column_name="NCBI_ID")
        genome_id_to_otu = metadata_table.get_map(key_column_name="genome_ID", value_column_name="OTU")

        column_genome_id = metadata_table.get_column("genome_ID")
        if not metadata_table.has_column("strain_id"):
            column_strain_id = metadata_table.get_empty_column()
        else:
            column_strain_id = metadata_table.get_column("strain_id")
            genome_id_to_strain_id = metadata_table.get_map(key_column_name="genome_ID", value_column_name="strain_id")

        genome_id_to_lineage = self._get_genome_id_to_lineage(
            genome_id_to_percent.keys(), genome_id_to_taxid, strain_id_to_genome_id, genome_id_to_strain_id)

        percent_by_rank_by_taxid = self._get_percent_by_rank_by_taxid(genome_id_to_lineage, genome_id_to_percent)

        # add strain_id to metadata
        #for row_index, genome_id in enumerate(column_genome_id):
        #    column_strain_id[row_index] = genome_id_to_strain_id[genome_id]
        #assert len(column_strain_id) == len(set(column_strain_id))
        #metadata_table.insert_column(column_strain_id, "strain_id")

        # stream taxonomic profile
        self._stream_tp_header(stream_output, sample_id)
        self._stream_tp_rows(stream_output, percent_by_rank_by_taxid, strain_id_to_genome_id, genome_id_to_otu)

    def _get_genome_id_to_lineage(
        self, list_of_genome_id, genome_id_to_taxid, strain_id_to_genome_id, genome_id_to_strain_id):
        """
        Returnes the lineage for each genome id, assigning new strain id if not available

        @param list_of_genome_id: List of identifier of genomes
        @type list_of_genome_id: list[str|unicode]
        @param genome_id_to_taxid: Assigned taxid for each genome id
        @type genome_id_to_taxid: dict[str|unicode, str|unicode]
        @param strain_id_to_genome_id: Mapping from strain id to genome id
        @type strain_id_to_genome_id: dict[str|unicode, str|unicode]
        @param genome_id_to_strain_id: Mapping from genome id to strain id
        @type genome_id_to_strain_id: dict[str|unicode, str|unicode]

        @return: lineage for each genome id using genome id as key
        @rtype: dict[str|unicode, list[None|str|unicode]]
        """
        strains_by_taxid = {}
        genome_id_to_lineage = {}
        for genome_id in list_of_genome_id:
            tax_id = genome_id_to_taxid[genome_id]
            if tax_id == "":
                raise KeyError("genome_ID '{}' has no taxid!".format(genome_id))
            tax_id = self._taxonomy.get_updated_taxid(tax_id)
            genome_id_to_lineage[genome_id] = self._taxonomy.get_lineage_of_legal_ranks(
                tax_id, ranks=self._ranks, default_value=None)
            if genome_id_to_lineage[genome_id][-1] is not None:
                continue

            if tax_id not in strains_by_taxid:
                strains_by_taxid[tax_id] = 0
            strains_by_taxid[tax_id] += 1

            if genome_id in genome_id_to_strain_id and genome_id_to_strain_id[genome_id]:
                strain_id = genome_id_to_strain_id[genome_id]
            else:
                strain_id = "{}.{}".format(tax_id, strains_by_taxid[tax_id])
                # make sure assigned strain ids are unique, in case of previous assigned ids
                while strain_id in genome_id_to_strain_id.values():
                    strains_by_taxid[tax_id] += 1
                    strain_id = "{}.{}".format(tax_id, strains_by_taxid[tax_id])
                genome_id_to_strain_id[genome_id] = strain_id
            genome_id_to_lineage[genome_id][-1] = strain_id
            strain_id_to_genome_id[strain_id] = genome_id
        return genome_id_to_lineage

    def _get_percent_by_rank_by_taxid(self, genome_id_to_lineage, genome_id_to_percent):
        """
        Return the percentage for each taxid of a list of default ranks

        @param genome_id_to_lineage: Mapping from genome id to a lineage (list)
        @type genome_id_to_lineage: dict[str|unicode, list[None|str|unicode]]
        @param genome_id_to_percent: Mapping from genome id to percentage
        @type genome_id_to_percent: dict[str|unicode, float]

        @return: Percentage for each taxid of a list of default ranks as dictionary of dictionaries
        @rtype: dict[str|unicode, dict[str|unicode, float]]
        """
        percent_by_rank_by_taxid = {}
        for rank in self._ranks:
            percent_by_rank_by_taxid[rank] = dict()

        for rank_index, rank in enumerate(self._ranks):
            # rank = ranks[rank_index]
            for genome_id in genome_id_to_lineage:
                tax_id = genome_id_to_lineage[genome_id][rank_index]
                if tax_id is None:
                    continue
                percent = genome_id_to_percent[genome_id]
                if tax_id not in percent_by_rank_by_taxid[rank]:
                    percent_by_rank_by_taxid[rank][tax_id] = 0
                percent_by_rank_by_taxid[rank][tax_id] += percent
        return percent_by_rank_by_taxid

    def _stream_tp_rows(self, stream_output, percent_by_rank_by_taxid, strain_id_to_genome_id, genome_id_to_otu):
        """
        Stream the rows of the taxonomic profile.

        @param stream_output: Output of taxonomic profile
        @type stream_output: file | FileIO | StringIO
        @param percent_by_rank_by_taxid: Percentage for each taxid of a list of default ranks as dictionary of dictionaries
        @type percent_by_rank_by_taxid: dict[str|unicode, dict[str|unicode, float]]
        @param strain_id_to_genome_id: Map from strain id to a genome identifier
        @type strain_id_to_genome_id: dict[str|unicode, str|unicode]
        @param genome_id_to_otu: Map from genome id to an otu identifier
        @type genome_id_to_otu: dict[str|unicode, str|unicode]
        """
        row_format = "{taxid}\t{rank}\t{taxpath}\t{taxpath_sn}\t{abp:.4f}\t{gid}\t{otu}\n"
        for rank_index, rank in enumerate(self._ranks):
            for tax_id in percent_by_rank_by_taxid[rank]:
                if tax_id == '':
                    self._logger.warning("Missing rank %s for a genome" % rank)
                    continue
                if '.' in tax_id:
                    genome_id = strain_id_to_genome_id[tax_id]
                    otu = genome_id_to_otu[genome_id]
                    lineage = self._taxonomy.get_lineage_of_legal_ranks(tax_id.split('.')[0], ranks=self._ranks, default_value="")
                    lineage[-1] = tax_id
                else:
                    genome_id = ""
                    otu = ""
                    lineage = self._taxonomy.get_lineage_of_legal_ranks(tax_id, ranks=self._ranks, default_value="")

                lineage = lineage[:rank_index+1]
                lineage_sn = [self._taxonomy.get_scientific_name(tid) if tid != "" and '.' not in tid else "" for tid in lineage]
                if '.' in tax_id:
                    lineage_sn[-1] = self._taxonomy.get_scientific_name(tax_id.split('.')[0]) + " strain"  # ""
                
                if percent_by_rank_by_taxid[rank][tax_id] != 0:
                    stream_output.write(row_format.format(
                        taxid=tax_id,
                        rank=rank,
                        taxpath="|".join(lineage),
                        taxpath_sn="|".join(lineage_sn),
                        abp=percent_by_rank_by_taxid[rank][tax_id]*100,
                        gid=genome_id,
                        otu=otu
                    ))

    def _stream_tp_header(self, output_stream, identifier):
        """
        Stream the header of the taxonomic profile.

        @param output_stream: Output of taxonomic profile
        @type output_stream: file | FileIO | StringIO
        @param identifier: Identifier of a sample
        @type identifier: str | unicode
        """
        output_stream.write("@SampleID:{}\n".format(identifier))
        output_stream.write("@Version:{}\n".format(self._taxonomic_profile_version))
        output_stream.write("@Ranks:{ranks}\n\n".format(ranks="|".join(self._ranks)))
        output_stream.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n")
