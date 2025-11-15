#!/usr/bin/env python

import sys
import os
import re

TAXID_TO_NAME = {}
TAXID_TO_PARENT_TAXID = {}
TAXID_TO_RANK = {}
TAXID_OLD_TO_TAXID_NEW = {}
LEGACY_RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
FILENAME_TAXONOMIC_PROFILE = "taxonomic_profile_{sample_index}.txt"
TAXONOMIC_PROFILE_VERSION = "0.9.2"

ROOTS = ["acellular root", "cellular root", "other entries"]
TOP   = ["realm", "domain"]
TAIL  = ["kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]

VIRUS_RANKS    = ['acellular root', 'realm',  'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
CELLULAR_RANKS = ['cellular root',  'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
OTHER_RANKS    = ['other entries',  'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']


# -----------------------------
# Constants for plasmids / other-entries
# -----------------------------
OTHER_ENTRIES_TAXID = "2787854"            # “other entries” synthetic root (if not in NCBI taxdump)
OTHER_ENTRIES_NAME  = "other entries"
PLASMID_SPECIES_TAXID = "45202"            # “unidentified plasmid” species
PLASMID_SPECIES_NAME  = "unidentified plasmid"

# read NCBI names file
def read_names_file(file_path_ncbi_names):
    with open(file_path_ncbi_names) as fin:
        for line in fin:
            # 65      |       Herpetosiphon aurantiacus       |               |       scientific name |
            taxid, name, _, nametype, *_ = line.strip().split('|')
            if nametype.strip() == 'scientific name':
                TAXID_TO_NAME[taxid.strip()] = name.strip()
    # inject names for our synthetic / special nodes if missing
    TAXID_TO_NAME.setdefault(OTHER_ENTRIES_TAXID, OTHER_ENTRIES_NAME)
    TAXID_TO_NAME.setdefault(PLASMID_SPECIES_TAXID, PLASMID_SPECIES_NAME)


def build_ncbi_taxonomy(file_path_ncbi_nodes):
    with open(file_path_ncbi_nodes) as f:
        for line in f:
            taxid, parent_taxid, rank, *_ = [el.strip() for el in line.split('|')]
            TAXID_TO_PARENT_TAXID[taxid] = parent_taxid
            TAXID_TO_RANK[taxid] = rank.lower()


# read NCBI merged file
def read_merged_file(file_path_ncbi_merged):
    with open(file_path_ncbi_merged) as fin:
        for line in fin:
            # 5085       |       746128  |
            old_taxid, new_taxid, *_ = line.strip().split('|')
            TAXID_OLD_TO_TAXID_NEW[old_taxid.strip()] = new_taxid.strip()


def write_taxonomic_profile_from_abundance_files(list_of_file_paths, metadata_file_path):
    """
    Write a taxonomic profile file for each relative abundance file

    @param list_of_file_paths: List of abundance file paths
    @type list_of_file_paths: list[str | unicode]
    """

    for file_path in list_of_file_paths:
        m = re.search(r'_([0-9]+)\.tsv$', os.path.basename(file_path))
        if not m:
            raise ValueError("Abundance file path does not end with _<number>.tsv: {}".format(file_path))
        index_abundance = int(m.group(1))
        community_abundance = parse_file(file_path)
        file_path_output = os.path.join("./", FILENAME_TAXONOMIC_PROFILE.format(sample_index=index_abundance))
        with open(file_path_output, 'w') as stream_output:
            write_taxonomic_profile(community_abundance, stream_output, index_abundance, metadata_file_path)


def write_taxonomic_profile(community_abundance, stream_output, sample_id, metadata_file_path):
    """
    Stream a taxonomic profile by list of relative abundances

    @param community_abundance: list of relative abundances
    @type community_abundance: generator[ dict[int|long|str|unicode, str|unicode] ]
    @param stream_output: Output of taxonomic profile
    @type stream_output: file | FileIO | StringIO
    @param sample_id: The sample ID of this taxonomy profile.
    @type output_stream: str
    """

    genome_abundance = {}
    total_abundance = 0.0
    for genome_id, abundance in community_abundance:
        if genome_id in genome_abundance:
            raise IOError("genome id '{}' is not unique!".format(genome_id))
        genome_abundance[genome_id] = float(abundance)
        total_abundance += genome_abundance[genome_id]

    for key, value in genome_abundance.items():
        genome_abundance[key] = (value / total_abundance) if total_abundance > 0 else 0.0

    stream_taxonomic_profile(stream_output, genome_abundance, sample_id, metadata_file_path)


def stream_taxonomic_profile(stream_output, genome_id_to_percent, sample_id, metadata_file_path):
    """
    Stream a taxonomic profile by list of percentages by genome id

    @param stream_output: Output of taxonomic profile
    @type stream_output: file | FileIO | StringIO
    @param genome_id_to_percent: Percentage for each genome id
    @type genome_id_to_percent: dict[str|unicode, float]
    @param sample_id: The sample ID of this taxonomy profile.
    @type output_stream: str
    """

    strain_id_to_genome_id = {}
    genome_id_to_strain_id = {}
    genome_id_to_ncbi_id = {}
    genome_id_to_otu = {}
    genome_id_is_plasmid = {}
    genome_id_host = {}

    # metadata: genome_ID \t OTU \t NCBI_ID \t novelty_category [\t strain_id (optional)]
    with open(metadata_file_path) as metadata:
        header = metadata.readline().rstrip('\n').split('\t')
        col = {h: i for i, h in enumerate(header)}

        # required columns
        required_cols = ["genome_ID", "OTU", "NCBI_ID", "novelty_category"]
        missing = [c for c in required_cols if c not in col]
        if missing:
            raise ValueError(
                "Metadata file is missing required columns: {}".format(", ".join(missing))
            )

        idx_genome = col["genome_ID"]
        idx_otu = col["OTU"]
        idx_ncbi = col["NCBI_ID"]
        idx_novelty = col["novelty_category"]
        idx_strain = col.get("strain_id", None)
        idx_host = col.get("host", None)
        max_idx = max(idx_genome, idx_otu, idx_ncbi, idx_novelty)

        for line in metadata:
            parts = line.rstrip('\n').split('\t')
            if not parts or len(parts) <= max_idx:
                continue

            genome_id = parts[idx_genome]
            strain_id = parts[idx_strain] if idx_strain is not None else None

            genome_id_to_ncbi_id[genome_id] = parts[idx_ncbi]
            genome_id_to_otu[genome_id] = parts[idx_otu]
            genome_id_is_plasmid[genome_id] = parts[idx_novelty].lower() == "plasmid"
            genome_id_host[genome_id] = parts[idx_host] if idx_host is not None else ''

            if strain_id:
                strain_id_to_genome_id[strain_id] = genome_id
                genome_id_to_strain_id[genome_id] = strain_id

    genome_id_to_lineage = get_genome_id_to_lineage(genome_id_to_percent.keys(), genome_id_to_ncbi_id,
                                                    strain_id_to_genome_id, genome_id_to_strain_id, genome_id_is_plasmid)

    percent_by_rank_by_taxid = get_percent_by_rank_by_taxid(genome_id_to_lineage, genome_id_to_percent, genome_id_is_plasmid)

    # stream taxonomic profile
    stream_tp_header(stream_output, sample_id)
    stream_tp_rows(stream_output, percent_by_rank_by_taxid, strain_id_to_genome_id, genome_id_to_otu, genome_id_host)


def stream_tp_rows(stream_output, percent_by_rank_by_taxid, strain_id_to_genome_id, genome_id_to_otu, genome_id_host):
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

    row_format = "{taxid}\t{rank}\t{taxpath}\t{taxpath_sn}\t{abp:.4f}\t{gid}\t{otu}\t{host}\n"
    all_rows = []

    ranks_map = {
        "acellular": VIRUS_RANKS,
        "cellular": CELLULAR_RANKS,
        "other": OTHER_RANKS
    }

    # Final output ordering
    master_rank_order = ROOTS + TOP + TAIL

    if "acellular root" not in set(TAXID_TO_RANK.values()):
        ranks_map = {
            "legacy": LEGACY_RANKS,
        }
        master_rank_order = LEGACY_RANKS

    # Collect rows
    for rank_type in ranks_map.keys():
        rank_list = ranks_map[rank_type]
        for rank in rank_list:
            taxid_to_percent = percent_by_rank_by_taxid.get(rank_type, {}).get(rank, {})
            for taxid, percent in taxid_to_percent.items():
                if not taxid:
                    continue

                is_strain = '.' in taxid
                if is_strain:
                    base_taxid = taxid.split('.')[0]
                    genome_id = strain_id_to_genome_id.get(taxid, "")
                    otu = genome_id_to_otu.get(genome_id, "")
                else:
                    base_taxid = taxid
                    genome_id = ""
                    otu = ""

                # Build lineage strings
                if rank_type == "other":  # plasmid track
                    lineage = get_other_entries_lineage(base_taxid, taxid_if_strain=taxid if is_strain else None)
                    lineage_sn = get_other_entries_lineage_names(base_taxid, taxid_if_strain=taxid if is_strain else None)
                else:
                    lineage = get_lineage_of_legal_ranks(base_taxid, default_value="")
                    lineage_sn = [
                        get_scientific_name(tid) if tid != "" and '.' not in tid else ""
                        for tid in lineage
                    ]
                    if is_strain:
                        lineage_sn[-1] = get_scientific_name(base_taxid) + " strain"
                        lineage[-1] = taxid

                host_str = genome_id_host.get(genome_id, "") if (rank_type == "other" and rank == "strain" and genome_id) else ""

                all_rows.append({
                    "taxid": taxid,
                    "rank": rank,
                    "rank_type": rank_type,
                    "percent": percent,
                    "rank_index": ranks_map[rank_type].index(rank),
                    "taxpath": lineage,
                    "taxpath_sn": lineage_sn,
                    "genome_id": genome_id,
                    "otu": otu,
                    "host": host_str
                })

    grouped_rows = {rank: [] for rank in master_rank_order}
    for row in all_rows:
        grouped_rows[row["rank"]].append(row)

    # Sort each group by abundance descending
    for rank in grouped_rows:
        grouped_rows[rank].sort(key=lambda r: -r["percent"])

    # Emit
    for rank in master_rank_order:
        for row in grouped_rows[rank]:
            stream_output.write(row_format.format(
                taxid=row["taxid"],
                rank=row["rank"],
                taxpath="|".join(row["taxpath"][:row["rank_index"] + 1]),
                taxpath_sn="|".join(row["taxpath_sn"][:row["rank_index"] + 1]),
                abp=row["percent"] * 100.0,
                gid=row["genome_id"],
                otu=row["otu"],
                host=row["host"]
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
    if taxid in TAXID_TO_NAME:
        return TAXID_TO_NAME[taxid]
    raise ValueError("Invalid taxid")


def stream_tp_header(output_stream, sample_id):
    """
    Stream the header of the taxonomic profile.

    @param output_stream: Output of taxonomic profile
    @type output_stream: file | FileIO | StringIO
    @param sample_id: The sample ID of this taxonomy profile.
    @type output_stream: str
    """
    output_stream.write("@SampleID:{}\n".format(sample_id))
    output_stream.write("@Version:{}\n".format(TAXONOMIC_PROFILE_VERSION))
    if "acellular root" not in set(TAXID_TO_RANK.values()):
        output_stream.write("@Ranks:{ranks}\n\n".format(ranks="|".join(LEGACY_RANKS)))
    else:
        output_stream.write("@Ranks:" + ",".join(ROOTS) + "|" + ",".join(TOP) + "|" + "|".join(TAIL) + "\n")
    output_stream.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\tHOST\n")


def get_percent_by_rank_by_taxid(genome_id_to_lineage, genome_id_to_percent, genome_id_is_plasmid):
    percent_by_rank_by_taxid = {
        "cellular":  {rank: {} for rank in CELLULAR_RANKS},
        "acellular": {rank: {} for rank in VIRUS_RANKS},
        "other":     {rank: {} for rank in OTHER_RANKS}
    }

    is_legacy = False
    if "acellular root" not in set(TAXID_TO_RANK.values()):
        is_legacy = True
        percent_by_rank_by_taxid = {
            "legacy": {rank: {} for rank in LEGACY_RANKS}
        }

    for genome_id, lineage in genome_id_to_lineage.items():
        if not lineage:
            continue
        base_taxid = lineage[0]  # first filled element of its rank list (root for that track)
        if base_taxid is None or base_taxid == "":
            continue

        # Decide the track
        if is_legacy:
            rank_type = "legacy"
            rank_list = LEGACY_RANKS
        elif genome_id_is_plasmid.get(genome_id, False):
            rank_type = "other"
            rank_list = OTHER_RANKS
        elif is_descendant_of(base_taxid, "10239"):  # viruses
            rank_type = "acellular"
            rank_list = VIRUS_RANKS
        else:
            rank_type = "cellular"
            rank_list = CELLULAR_RANKS

        percent = genome_id_to_percent[genome_id]

        for i, taxid_i in enumerate(lineage):
            if taxid_i is None or taxid_i == "":
                continue
            rank = rank_list[i]
            if taxid_i not in percent_by_rank_by_taxid[rank_type][rank]:
                percent_by_rank_by_taxid[rank_type][rank][taxid_i] = 0.0
            percent_by_rank_by_taxid[rank_type][rank][taxid_i] += percent

    return percent_by_rank_by_taxid


def get_genome_id_to_lineage(list_of_genome_id, genome_id_to_taxid, strain_id_to_genome_id,
                             genome_id_to_strain_id, genome_id_is_plasmid):
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
        if genome_id_is_plasmid.get(genome_id, False):
            # Build OTHER-ENTRIES lineage
            lineage = ["" for _ in OTHER_RANKS]
            lineage[0] = OTHER_ENTRIES_TAXID
            lineage[8] = PLASMID_SPECIES_TAXID  # species slot
            # assign strain
            if genome_id in genome_id_to_strain_id:
                lineage[9] = genome_id_to_strain_id[genome_id]
            else:
                # synthesize a unique 45202.N strain id
                n = strains_by_taxid.get(PLASMID_SPECIES_TAXID, 0) + 1
                strains_by_taxid[PLASMID_SPECIES_TAXID] = n
                sid = f"{PLASMID_SPECIES_TAXID}.{n}"
                genome_id_to_strain_id[genome_id] = sid
                strain_id_to_genome_id[sid] = genome_id
                lineage[9] = sid
            genome_id_to_lineage[genome_id] = lineage
            continue

        # regular (cellular or viral) genomes
        tax_id = genome_id_to_taxid[genome_id]
        if tax_id == "":
            raise KeyError("genome_ID '{}' has no taxid!".format(genome_id))
        tax_id = get_updated_taxid(tax_id)

        lineage = get_lineage_of_legal_ranks(tax_id)
        # Ensure a strain id exists
        if lineage[-1] is None:
            if tax_id not in strains_by_taxid:
                strains_by_taxid[tax_id] = 0
            strains_by_taxid[tax_id] += 1
            strain_id = genome_id_to_strain_id.get(genome_id, "{}.{}".format(tax_id, strains_by_taxid[tax_id]))
            while strain_id in strain_id_to_genome_id:
                strains_by_taxid[tax_id] += 1
                strain_id = "{}.{}".format(tax_id, strains_by_taxid[tax_id])
            genome_id_to_strain_id[genome_id] = strain_id
            strain_id_to_genome_id[strain_id] = genome_id
            lineage[-1] = strain_id

        genome_id_to_lineage[genome_id] = lineage
    return genome_id_to_lineage


def get_lineage_of_legal_ranks(taxid, default_value=None):
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

    taxid = get_updated_taxid(taxid)
    is_virus = is_descendant_of(taxid, "10239")
    ranks = VIRUS_RANKS if is_virus else CELLULAR_RANKS
    lineage_dict = {rank: default_value for rank in ranks}
    while taxid != "1":
        if taxid not in TAXID_TO_RANK or taxid not in TAXID_TO_PARENT_TAXID:
            break
        rank = TAXID_TO_RANK[taxid]
        if rank in lineage_dict:
            lineage_dict[rank] = taxid
        taxid = TAXID_TO_PARENT_TAXID[taxid]
    return [lineage_dict[rank] for rank in ranks]


def get_updated_taxid(taxid):
    """
    Return current taxid, in case it was merged

    @attention: taxid is not accepted as digit!!!

    @param taxid: ncbi taxonomic identifier
    @type taxid: str

    @return: ncbi taxonomic identifier
    @rtype: str | unicode
    """

    if taxid in TAXID_TO_RANK:
        return taxid
    if taxid in TAXID_OLD_TO_TAXID_NEW:
        return TAXID_OLD_TO_TAXID_NEW[taxid]
    # allow our injected nodes
    if taxid in (OTHER_ENTRIES_TAXID, PLASMID_SPECIES_TAXID):
        return taxid
    raise ValueError(f"Invalid taxid: {taxid}")


def parse_file(file_path):
    """
	Reading comma or tab separated values from a file

	@param file_path: path to file to be opened
	@type file_path: str | unicode

	@return: Generator of dictionary representing rows
	@rtype: generator[ dict[int|long|str|unicode, str|unicode] ]
	"""

    with open(file_path) as file_handler:
        for row in parse_stream(file_handler):
            yield row


def parse_stream(stream_input):
    """
	Reading comma or tab separated values from a stream

	@param stream_input: stream
	@type stream_input: file | io.FileIO | StringIO.StringIO
    """

    separator = "\t"
    for line in stream_input:
        row = line.rstrip('\n').rstrip('\r')
        if not row or line[0] == '#':
            continue
        row_cells = row.split(separator)
        yield row_cells


def get_other_entries_lineage(species_taxid, taxid_if_strain=None):
    """Build lineage aligned to other_ranks, with OTHER_ENTRIES_TAXID at root and
    species_taxid at 'species'. Optionally set 'strain' to taxid_if_strain."""
    lineage = [""] * len(OTHER_RANKS)
    lineage[0] = OTHER_ENTRIES_TAXID
    lineage[8] = species_taxid
    if taxid_if_strain:
        lineage[9] = taxid_if_strain
    return lineage


def get_other_entries_lineage_names(species_taxid, taxid_if_strain=None):
    names = ["" for _ in OTHER_RANKS]
    names[0] = OTHER_ENTRIES_NAME
    names[8] = PLASMID_SPECIES_NAME if species_taxid == PLASMID_SPECIES_TAXID else TAXID_TO_NAME.get(species_taxid, "")
    if taxid_if_strain:
        names[9] = (PLASMID_SPECIES_NAME + " strain")
    return names


def is_descendant_of(taxid, ancestor_taxid):
    seen = set()
    while taxid != "1" and taxid not in seen:
        seen.add(taxid)
        if taxid == ancestor_taxid:
            return True
        if taxid not in TAXID_TO_PARENT_TAXID:
            break
        taxid = TAXID_TO_PARENT_TAXID[taxid]
    return False


# main method and entry point of this script
# this script builds the taxonomy for a given ncbi dump and given distributions
if __name__ == "__main__":
    file_path_ncbi_names = sys.argv[1]
    file_path_ncbi_merged = sys.argv[2]
    file_path_ncbi_nodes = sys.argv[3]
    sample_size = int(sys.argv[4])
    metadata_file_path = sys.argv[5]
    list_of_file_paths_distribution = [sys.argv[6 + i] for i in range(sample_size)]

    build_ncbi_taxonomy(file_path_ncbi_nodes)
    read_names_file(file_path_ncbi_names)
    read_merged_file(file_path_ncbi_merged)

    write_taxonomic_profile_from_abundance_files(list_of_file_paths_distribution, metadata_file_path)