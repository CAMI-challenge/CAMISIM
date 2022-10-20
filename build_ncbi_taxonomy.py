#!/usr/bin/env python

import sys 
import os
import shutil

taxid_to_name = {}
taxid_to_parent_taxid = {}
taxid_to_rank = {}
taxid_old_to_taxid_new = {}
ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
filename_taxonomic_profile = "taxonomic_profile_{sample_index}.txt"
default_ordered_legal_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
metadata_file_path = "/home/jfunk/CAMISIM/code/CAMISIM/defaults/metadata.tsv"
taxonomic_profile_version = "0.9.1"

# read NCBI names file
def read_names_file(file_path_ncbi_names):

    with open(file_path_ncbi_names) as fin:
            for line in fin:
                # 65      |       Herpetosiphon aurantiacus       |               |       scientific name |
                taxid, name, disambiguation, nametype, more = line.strip().split('|')
                if nametype.strip() == 'scientific name':
                    taxid_to_name[taxid.strip()] = name.strip()


def build_ncbi_taxonomy(file_path_ncbi_nodes):

    with open(file_path_ncbi_nodes) as file_handler:
            for line in file_handler:
                elements = [el.strip() for el in line.split('|')]
                taxid, parent_taxid, rank = elements[0:3]
                rank = rank.lower()  # should be lower-case in file, but can't be bad to doublecheck
                taxid_to_parent_taxid[taxid] = parent_taxid
                taxid_to_rank[taxid] = rank

# read NCBI merged file
def read_merged_file(file_path_ncbi_merged):

    with open(file_path_ncbi_merged) as fin:
            for line in fin:
                # 5085       |       746128  |
                old_taxid, new_taxid, sonst = line.strip().split('|')
                taxid_old_to_taxid_new[old_taxid.strip()] = new_taxid.strip()

def write_taxonomic_profile_from_abundance_files(list_of_file_paths):
    
    """
    Write a taxonomic profile file for each relative abundance file

    @param list_of_file_paths: List of abundance file paths
    @type list_of_file_paths: list[str | unicode]
    """
    
    for index_abundance, file_path in enumerate(list_of_file_paths):
        community_abundance = parse_file(file_path)
        file_path_output = os.path.join("./", filename_taxonomic_profile.format(
            sample_index=index_abundance))
        with open(file_path_output, 'w') as stream_output:
            write_taxonomic_profile(community_abundance, stream_output)


def write_taxonomic_profile(community_abundance, stream_output):

    """
    Stream a taxonomic profile by list of relative abundances

    @param community_abundance: list of relative abundances
    @type community_abundance: generator[ dict[int|long|str|unicode, str|unicode] ]
    @param stream_output: Output of taxonomic profile
    @type stream_output: file | FileIO | StringIO
    """

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

    stream_taxonomic_profile(stream_output, genome_abundance)


def stream_taxonomic_profile(stream_output, genome_id_to_percent):

    """
    Stream a taxonomic profile by list of percentages by genome id

    @param stream_output: Output of taxonomic profile
    @type stream_output: file | FileIO | StringIO
    @param genome_id_to_percent: Percentage for each genome id
    @type genome_id_to_percent: dict[str|unicode, float]
    """

    strain_id_to_genome_id = {}
    genome_id_to_strain_id = {}
    genome_id_to_ncbi_id = {}
    genome_id_to_otu = {}

    counter = 0
    strain_id_in_metadata = False

    with open(metadata_file_path) as metadata:
        for line in metadata:
            if(counter == 0):
                column_names = line.strip().split('\t')
                if("strain_id" in column_names):
                    strain_id_in_metadata = True
            else:
                if(strain_id_in_metadata):
                    genome_id, otu, ncbi_id, novelty_category, strain_id = line.strip().split('\t')
                    strain_id_to_genome_id[strain_id] = genome_id
                    genome_id_to_strain_id[genome_id] = strain_id
                else:
                    genome_id, otu, ncbi_id, novelty_category = line.strip().split('\t')
                genome_id_to_ncbi_id[genome_id] = ncbi_id
                genome_id_to_otu[genome_id] = otu
            counter = counter + 1

    genome_id_to_lineage = get_genome_id_to_lineage(genome_id_to_percent.keys(), genome_id_to_ncbi_id, strain_id_to_genome_id, genome_id_to_strain_id)

    percent_by_rank_by_taxid = get_percent_by_rank_by_taxid(genome_id_to_lineage, genome_id_to_percent)

    # add strain_id to metadata
    #for row_index, genome_id in enumerate(column_genome_id):
    #    column_strain_id[row_index] = genome_id_to_strain_id[genome_id]
    #assert len(column_strain_id) == len(set(column_strain_id))
    #metadata_table.insert_column(column_strain_id, "strain_id")

    # stream taxonomic profile
    stream_tp_header(stream_output)
    stream_tp_rows(stream_output, percent_by_rank_by_taxid, strain_id_to_genome_id, genome_id_to_otu)


def stream_tp_rows(stream_output, percent_by_rank_by_taxid, strain_id_to_genome_id, genome_id_to_otu):

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
    for rank_index, rank in enumerate(ranks):
        for tax_id in percent_by_rank_by_taxid[rank]:
            if tax_id == '':
                continue
            if '.' in tax_id:
                genome_id = strain_id_to_genome_id[tax_id]
                otu = genome_id_to_otu[genome_id]
                lineage = get_lineage_of_legal_ranks(tax_id.split('.')[0], ranks=ranks, default_value="")
                lineage[-1] = tax_id
            else:
                genome_id = ""
                otu = ""
                lineage = get_lineage_of_legal_ranks(tax_id, ranks=ranks, default_value="")

            lineage = lineage[:rank_index+1]
            lineage_sn = [get_scientific_name(tid) if tid != "" and '.' not in tid else "" for tid in lineage]
            if '.' in tax_id:
                lineage_sn[-1] = get_scientific_name(tax_id.split('.')[0]) + " strain"  # ""
                
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


def get_scientific_name(taxid):

    """
    Return scientific name of ncbi taxonomic identifier

    @attention: taxid is not accepted as digit!!!

    @param taxid: ncbi taxonomic identifier
    @type taxid: str

    @return: ncbi scientific name
    @rtype: str | unicode
    """

    taxid = get_updated_taxid(taxid)
    if taxid in taxid_to_name:
        return taxid_to_name[taxid]
    raise ValueError("Invalid taxid")


def stream_tp_header(output_stream):
        """
        Stream the header of the taxonomic profile.

        @param output_stream: Output of taxonomic profile
        @type output_stream: file | FileIO | StringIO
        """
        identifier = ""
        output_stream.write("@SampleID:{}\n".format(identifier))
        output_stream.write("@Version:{}\n".format(taxonomic_profile_version))
        output_stream.write("@Ranks:{ranks}\n\n".format(ranks="|".join(ranks)))
        output_stream.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n")


def get_percent_by_rank_by_taxid(genome_id_to_lineage, genome_id_to_percent):

    percent_by_rank_by_taxid = {}
    for rank in ranks:
        percent_by_rank_by_taxid[rank] = dict()

    for rank_index, rank in enumerate(ranks):
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


def get_genome_id_to_lineage(list_of_genome_id, genome_id_to_taxid, strain_id_to_genome_id, genome_id_to_strain_id):

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
        tax_id = get_updated_taxid(tax_id)
        genome_id_to_lineage[genome_id] = get_lineage_of_legal_ranks(tax_id, default_ordered_legal_ranks)
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


def get_lineage_of_legal_ranks(taxid, ranks, default_value=None, as_name=False, inherit_rank=False):

    """
    Return lineage of a specific taxonomic identifier, filtered by a list of legal ranks

    @attention: The list of ranks determines the order of the returned taxonomic identifiers

    @param taxid: ncbi taxonomic identifier
    @type taxid: str
    @param ranks: List of ncbi ranks in lower case
    @type ranks: list[str]
    @param default_value: Value at rank indexes at which the taxid of that specific rank is undefined
    @type default_value: None | str
    @param as_name: return scientific name if true, not taxonomic id
    @type as_name: bool
    @param inherit_rank: name unnamed rank names by known ones, species -> root
    @type inherit_rank: bool

    @return: list of ncbi taxonomic identifiers
    @rtype: list[str|unicode|None]
    """

    assert isinstance(taxid, str)
    taxid = get_updated_taxid(taxid)

    lineage = [default_value] * len(ranks)
    original_rank = get_rank_of_taxid(taxid)
    if original_rank is not None and original_rank in ranks:
        if as_name:
            lineage[ranks.index(original_rank)] = taxid_to_name[taxid]
        else:
            lineage[ranks.index(original_rank)] = taxid
    try:
        rank_counter = ranks.index(taxid_to_rank[taxid]) # starting at rank of original tax id
    except ValueError: # rank is not in ranks
        rank_counter = ranks.index(ranks[-1]) # choose lowest rank then
    while taxid != "1":
        taxid = taxid_to_parent_taxid[taxid]
        rank = taxid_to_rank[taxid]
        if rank in ranks:
            current_rank_counter = ranks.index(taxid_to_rank[taxid])
            rank_difference = rank_counter - current_rank_counter
            if rank_difference > 1:
                for i in range(current_rank_counter, rank_counter - 1):
                    lineage[i] = "" # add empty name to list if name is missing in the taxonomy
            rank_counter = current_rank_counter
            if as_name:
                lineage[ranks.index(rank)] = taxid_to_name[taxid]
            else:
                lineage[ranks.index(rank)] = taxid

        # todo: sort ranks
        if inherit_rank:
            rank_previous = default_value
            tmp_list = enumerate(lineage)
            if default_ordered_legal_ranks.index(ranks[0]) < default_ordered_legal_ranks.index(ranks[-1]):
                tmp_list = reversed(list(enumerate(lineage)))
            for index, value in tmp_list:
                if value == default_value:
                    lineage[index] = rank_previous
                else:
                    rank_previous = value
    return lineage

def get_rank_of_taxid(taxid):

    """
    Return rank of ncbi taxonomic identifier

    @param taxid: ncbi taxonomic identifier
    @type taxid: str

    @return: ncbi rank of taxonomic identifiers
    @rtype: str | unicode
    """

    assert isinstance(taxid, str)
    taxid = get_updated_taxid(taxid)
    if taxid in taxid_to_rank:
        return taxid_to_rank[taxid]
    raise ValueError("Invalid taxid")

def get_updated_taxid(taxid):

    """
    Return current taxid, in case it was merged

    @attention: taxid is not accepted as digit!!!

    @param taxid: ncbi taxonomic identifier
    @type taxid: str

    @return: ncbi taxonomic identifier
    @rtype: str | unicode
    """

    if taxid in taxid_to_rank:
        return taxid
    if taxid not in taxid_old_to_taxid_new:
        raise ValueError("Invalid taxid")

    taxid_new = taxid_old_to_taxid_new[taxid]
    return taxid_new


def parse_file(file_path):

    """
	Reading comma or tab separated values from a file

	@param file_path: path to file to be opened
	@type file_path: str | unicode

	@return: Generator of dictionary representing rows
	@rtype: generator[ dict[int|long|str|unicode, str|unicode] ]
	"""

    print (file_path)
    with open(file_path) as file_handler:
        for row in parse_stream(file_handler):
            yield row


def parse_stream(stream_input):

    """
	Reading comma or tab separated values from a stream

	@param stream_input: stream
	@type stream_input: file | io.FileIO | StringIO.StringIO
    """

    comment_line = ['#']
    separator="\t"

    # read column names
    number_of_columns = 0
    list_of_column_names = []
	
    # read rows
    line_count = 0
    for line in stream_input:
        line_count += 1
        row = line.rstrip('\n').rstrip('\r')
        if line[0] in comment_line or len(row) == 0:
            continue

        row_cells = row.split(separator)
        if number_of_columns == 0:
            number_of_columns = len(row_cells)

        if number_of_columns != len(row_cells):
            msg = "Format error. Bad number of values in line {}".format(line_count)
            raise ValueError(msg)

        yield row_cells

# main method and entry point of this script
# this script builds the taxonomy for a given ncbi dump and given distributions
if __name__ == "__main__":
    file_path_ncbi_names = sys.argv[1]
    file_path_ncbi_merged = sys.argv[2]
    file_path_ncbi_nodes = sys.argv[3]
    sample_size = int(sys.argv[4])
    list_of_file_paths_distribution = list()
    i = 0

    for i in range(sample_size):
        list_of_file_paths_distribution.append(sys.argv[5+i])
        i = i + 1 


    build_ncbi_taxonomy(file_path_ncbi_nodes)
    read_names_file(file_path_ncbi_names)
    read_merged_file(file_path_ncbi_merged)

    write_taxonomic_profile_from_abundance_files(list_of_file_paths_distribution)