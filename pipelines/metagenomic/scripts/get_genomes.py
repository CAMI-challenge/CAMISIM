import sys
import biom
import os
import gzip
import urllib.request
import shutil
from numpy import random as np_rand
from ete4 import NCBITaxa
import sqlite3

number_of_samples = 0
LEGACY_RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
NEW_RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'realm', 'domain', 'acellular root', 'cellular root']
RANKS = []


"""
Reads a BIOM file and creates map of OTU: lineage, abundance
BIOM file format needs to have a taxonomy field in metadata which contains the taxonomy in the format:
RANK__SCINAME; LOWERRANK_LOWERSCINAME
"""
def read_taxonomic_profile(biom_profile, no_samples = None):
    table = biom.load_table(biom_profile)
    ids = table.ids(axis="observation")
    samples = table.ids()

    if no_samples is None:
        no_samples = len(samples)

    if no_samples is not None and no_samples != len(samples) and no_samples != 1:
        # log.warning("Number of samples (%s) does not match number of samples in biom file (%s)" % (no_samples, len(samples)))
        if no_samples > len(samples):
            no_samples = len(samples)
        #_log.warning("Using the first %s samples" % no_samples)

    profile = {}
    for otu in ids:
        lineage = table.metadata(otu,axis="observation")["taxonomy"]
        try:
            lineage = lineage.split(";") # if no spaces
        except AttributeError:
            pass
        abundances = []
        for sample in samples[:no_samples]:
            abundances.append(table.get_value_by_ids(otu,sample))
        profile[otu] = (lineage, abundances)
    
    return profile


def normalize_genome_key(genome_path):
    genome_path = genome_path.strip()
    if "://" in genome_path:
        return genome_path.rstrip("/")
    return os.path.normpath(genome_path)


def get_required_quality_column(columns, aliases):
    for alias in aliases:
        if alias in columns:
            return columns[alias]
    raise ValueError("Required column(s) missing from additional genome quality file: %s" % ", ".join(aliases))


def read_additional_genomes_quality(quality_file):
    quality_by_path = {}
    with open(quality_file, 'r') as quality_stream:
        header = quality_stream.readline().rstrip('\n').split('\t')
        if not header or header == ['']:
            raise ValueError("Additional genome quality file '%s' is empty" % quality_file)

        columns = {column.strip().lower(): idx for idx, column in enumerate(header)}
        idx_path = get_required_quality_column(columns, ['genome_path', 'path', 'genome'])
        idx_completeness = get_required_quality_column(columns, ['completeness'])
        idx_contamination = get_required_quality_column(columns, ['contamination'])
        idx_num_contigs = get_required_quality_column(columns, ['num_contigs', 'contigs'])
        max_idx = max(idx_path, idx_completeness, idx_contamination, idx_num_contigs)

        for line in quality_stream:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) <= max_idx:
                raise ValueError("Malformed row in additional genome quality file '%s': %s" % (quality_file, line))

            genome_key = normalize_genome_key(parts[idx_path])
            if genome_key in quality_by_path:
                raise ValueError("Duplicate quality entry for additional genome '%s'" % parts[idx_path])

            quality_by_path[genome_key] = {
                "completeness": float(parts[idx_completeness]),
                "contamination": float(parts[idx_contamination]),
                "num_contigs": int(float(parts[idx_num_contigs]))
            }
    return quality_by_path


def meets_additional_quality_requirements(quality_entry, max_contamination, min_completeness, max_num_contigs):
    return (
        quality_entry["contamination"] <= max_contamination and
        quality_entry["completeness"] >= min_completeness and
        quality_entry["num_contigs"] <= max_num_contigs
    )


"""
Reads list of available genomes in the (tsv) format:
NCBI_ID Scientific_Name ftp_path
Additional files might be provided with:
NCBI_ID Scientific_Name genome_path novelty_category
were path might either be online or offline/local
"""
def read_genomes_list(genomes_path, additional_file = None, additional_quality_file = None,
                      max_contamination = 5, min_completeness = 95, max_num_contigs = 200):
    genomes_map = {}
    total_genomes = 0
    additional_quality = {}
    if additional_file is not None:
        if additional_quality_file is None:
            raise ValueError("Parameter additional_genomes_quality_file is required when additional_references is set")
        additional_quality = read_additional_genomes_quality(additional_quality_file)
        with open(additional_file,'r') as add:
            for line in add:
                if not line.strip():
                    continue
                ncbi_id, sci_name, path, novelty = line.strip().split('\t')
                path_key = normalize_genome_key(path)
                if path_key not in additional_quality:
                    raise ValueError("Missing quality entry for additional genome '%s' in '%s'" % (path, additional_quality_file))
                quality_entry = additional_quality[path_key]
                genome_record = {
                    "path": path,
                    "novelty": novelty,
                    "source": "additional",
                    "quality_pass": meets_additional_quality_requirements(
                        quality_entry,
                        max_contamination,
                        min_completeness,
                        max_num_contigs
                    )
                }
                if ncbi_id in genomes_map:
                    genomes_map[ncbi_id]["genomes"].append(genome_record)
                else:
                    genomes_map[ncbi_id] = {"scientific_name": sci_name, "genomes": [genome_record]}
                total_genomes += 1
    with open(genomes_path,'r') as genomes:
        for line in genomes:
            if not line.strip():
                continue
            ncbi_id, sci_name, ftp = line.strip().split('\t')
            http = ftp.replace("ftp://","http://") # not using ftp address but http (proxies)
            genome_record = {
                "path": http,
                "novelty": 'known_strain',
                "source": "reference",
                "quality_pass": True
            }
            if ncbi_id in genomes_map:
                genomes_map[ncbi_id]["genomes"].append(genome_record)
            else:
                genomes_map[ncbi_id] = {"scientific_name": sci_name, "genomes": [genome_record]} # sci_name is always the same for same taxid (?)
            total_genomes += 1
    return genomes_map, total_genomes

"""
Given all available genomes, creates a map sorted by ranks of available genomes on that particular rank, ordered by their ncbi ids
"""
def get_genomes_per_rank(genomes_map, ncbi):
    per_rank_map = {}
    for rank in RANKS:
        per_rank_map[rank] = {}
    for genome in genomes_map:
        try:
            lineage = ncbi.get_lineage(genome) # this might contain some others ranks than ranks
            ranks_lin = ncbi.get_rank(lineage)
            genome_records = []
            for genome_record in genomes_map[genome]["genomes"]:
                genome_records.append({
                    "path": genome_record["path"],
                    "genome_id": genome,
                    "novelty": genome_record["novelty"],
                    "source": genome_record["source"],
                    "quality_pass": genome_record["quality_pass"]
                })
            for tax_id in lineage: # go over the lineage
                try:
                    check_rank = ranks_lin[tax_id]
                except KeyError:
                    continue
                if check_rank in per_rank_map: # if we are a legal rank
                    rank_map = per_rank_map[ranks_lin[tax_id]]
                    if tax_id not in rank_map: # tax id doesn't have a genome yet
                        rank_map[tax_id] = []
                    for genome_record in genome_records:
                        rank_map[tax_id].append(genome_record)
        except ValueError as e:
#           log.warning(e)
            print(e) # ToDo
    return per_rank_map

"""
Sorts the otus in the profile by abundance
"""
def sort_by_abundance(profile):
    sorted_keys = []
    for otu in profile:
        lineage, abundances = profile[otu]
        avg_abundance = sum(abundances)/len(abundances)
        sorted_keys.append((avg_abundance, otu)) # average abundance
    sorted_keys = sorted(sorted_keys, reverse=True)
    return [key for ab,key in sorted_keys]


"""
Given a BIOM lineage, create a NCBI tax id lineage
"""
def transform_lineage(lineage, ncbi):
    new_lineage = []
    for member in lineage:
        name = member.split("__")[-1] # name is on the right hand side
        if len(name) == 0:
            continue
        mapping = ncbi.get_name_translator([name])
        if name in mapping:
            taxid = mapping[name][0]
            if ncbi.get_rank([taxid])[taxid] in RANKS:
                new_lineage.append(taxid) # should contain only one element
        else:
            fallback_name = name.split()[0]
            fallback_mapping = ncbi.get_name_translator([fallback_name])
            if fallback_name in fallback_mapping:
                taxid = fallback_mapping[fallback_name][0]
                if ncbi.get_rank([taxid])[taxid] in RANKS:
                    new_lineage.append(taxid) # retry if space in name destroys ID
    return new_lineage[::-1] # invert list, so lowest rank appears first (last in BIOM)


def truncated_geometric(p, l, u, max_tries=100000):
    if not (0 < p <= 1) or l > u:
        raise ValueError("Invalid parameters: p must be in (0,1) and l <= u.")

    for _ in range(max_tries):
        k = np_rand.geometric(p)
        if l <= k <= u:
            return k

    raise RuntimeError(
        f"Failed to draw a sample in the range [{l}, {u}] within {max_tries} tries. "
        f"The probability of a draw falling in this range may be too low for p={p}."
    )


def set_abundances(otu_genome_map, abundances, used_genomes, otu, tax_id, mu, sigma):
    log_normal_vals = np_rand.lognormal(mu, sigma, len(used_genomes))
    sum_log_normal = sum(log_normal_vals)
    for i, genome in enumerate(used_genomes):
        otu_id = otu + "." + str(i)
        otu_genome_map[otu_id] = (tax_id, genome["genome_id"], genome["path"], genome["novelty"], [])  # taxid, genomeid, http path, novelty, abundances per sample
        relative_abundance = log_normal_vals[i] / sum_log_normal
        for abundance in abundances:  # calculate abundance per sample
            current_abundance = relative_abundance * abundance
            otu_genome_map[otu_id][-1].append(current_abundance)


def sample_genomes(genomes, amount):
    if amount <= 0 or not genomes:
        return []
    draw_count = min(len(genomes), amount)
    used_indices = np_rand.choice(len(genomes), draw_count, replace=False)
    return [genomes[i] for i in used_indices]


def get_candidate_groups(available_genomes, prioritize_additional_genomes):
    high_quality_additional = [genome for genome in available_genomes if genome["source"] == "additional" and genome["quality_pass"]]
    reference_genomes = [genome for genome in available_genomes if genome["source"] == "reference"]
    low_quality_additional = [genome for genome in available_genomes if genome["source"] == "additional" and not genome["quality_pass"]]

    if prioritize_additional_genomes:
        candidate_groups = [high_quality_additional, reference_genomes]
    else:
        candidate_groups = [high_quality_additional + reference_genomes]

    candidate_groups.append(low_quality_additional)
    return candidate_groups


def pick_genomes(tax_id, max_rank, min_strains, max_strains, per_rank_map, ncbi, no_replace, prioritize_additional_genomes):
    rank = ncbi.get_rank([tax_id])[tax_id]
    if RANKS.index(rank) > RANKS.index(max_rank):
        return []
    available_genomes = per_rank_map[rank].get(tax_id)
    if not available_genomes:
        return []
    strains_to_draw = truncated_geometric(2. / max_strains, min_strains, max_strains)
    used_genomes = []
    for candidate_group in get_candidate_groups(available_genomes, prioritize_additional_genomes):
        remaining = strains_to_draw - len(used_genomes)
        if remaining <= 0:
            break
        used_genomes.extend(sample_genomes(candidate_group, remaining))

    if no_replace: # sampling without replacement:
        for genome in used_genomes:
            for new_rank in per_rank_map:
                for taxid in per_rank_map[new_rank]:
                    if genome in per_rank_map[new_rank][taxid]:
                        per_rank_map[new_rank][taxid].remove(genome)

    return used_genomes


"""
Given the OTU to lineage/abundances map and the genomes to lineage map, create map otu: taxid, genome, abundances
"""
def map_otus_to_genomes(tax_profile, per_rank_map, mu, sigma, min_strains, max_strains, no_replace,
                        max_genomes, max_rank, prioritize_additional_genomes, ncbi):
    unmatched_otus = []
    matched_otus = set()
    otu_genome_map = {}
    sorted_otus = sort_by_abundance(tax_profile)
    genome_set_size = 0

    otu_to_lineage = {}
    otus_list = []
    max_lineage_len = 0
    for otu in sorted_otus:
        lin, _ = tax_profile[otu]
        otu_to_lineage[otu] = transform_lineage(lin, ncbi)
        if len(otu_to_lineage[otu]) == 0:
            unmatched_otus.append(otu)
            continue
        otus_list.append(otu)
        max_lineage_len = max(max_lineage_len, len(otu_to_lineage[otu]))

    for idx_lineage in range(max_lineage_len + 1):
        for otu in otus_list:
            if genome_set_size >= max_genomes and no_replace:  # cancel if no genomes are available anymore
                break
            if otu in matched_otus:
                continue

            try:
                tax_id = otu_to_lineage[otu][idx_lineage]
            except IndexError:
                continue

            used_genomes = pick_genomes(tax_id, max_rank, min_strains, max_strains, per_rank_map, ncbi, no_replace,
                                        prioritize_additional_genomes)
            if not used_genomes:
                continue

            genome_set_size += len(used_genomes)  # how many genomes are used

            matched_otus.add(otu)
            set_abundances(otu_genome_map, tax_profile[otu][1], used_genomes, otu, tax_id, mu, sigma)

    for otu in otus_list:
        if otu not in matched_otus:
            unmatched_otus.append(otu)

    return otu_genome_map, unmatched_otus

def fill_up_genomes(otu_genome_map, unmatched_otus, per_rank_map, tax_profile, debug, prioritize_additional_genomes):
    genomes = {}
    for rank in per_rank_map:
        for taxid in per_rank_map[rank]:
            for genome in per_rank_map[rank][taxid]:
                genome_key = (genome["path"], genome["genome_id"])
                if genome_key not in genomes:
                    genomes[genome_key] = {
                        "tax_id": taxid,
                        "genome_id": genome["genome_id"],
                        "path": genome["path"],
                        "novelty": genome["novelty"],
                        "source": genome["source"],
                        "quality_pass": genome["quality_pass"]
                    }

    ordered_genomes = []
    for candidate_group in get_candidate_groups(list(genomes.values()), prioritize_additional_genomes):
        ordered_genomes.extend(sample_genomes(candidate_group, len(candidate_group)))

    if not ordered_genomes:
        return otu_genome_map

    otu_indices = np_rand.choice(len(unmatched_otus),len(unmatched_otus),replace=False)
    for i, genome in enumerate(ordered_genomes[:len(unmatched_otus)]):
        curr_otu = unmatched_otus[otu_indices[i]] #so we choose a random genome
        lineage, abundances = tax_profile[curr_otu]
        otu_genome_map[curr_otu] = (genome["tax_id"], genome["genome_id"], genome["path"], genome["novelty"], abundances)
#        if debug:
#            _log.warning("Filling up OTU %s (mapped tax id: %s) to genome with tax id %s" % (curr_otu, lin[0], tax_id))
    return otu_genome_map

"""
Take fasta input file and split by any N occurence (and remove Ns)
"""
def split_by_N(fasta_path, out_path, script):
    #os.system("scripts/split_fasta.pl %s %s" % (fasta_path, out_path))
    os.system(script+" %s %s" % (fasta_path, out_path))
    os.remove(fasta_path)

"""
Downloads the given genome and returns the out path
"""
def download_genome(genome, out_path, script):
    # genome_path = os.path.join(out_path,"genomes")
    genome_path = out_path


    out_name = genome.rstrip().split('/')[-1]
    http_address = os.path.join(genome, out_name + "_genomic.fna.gz")
    opened = urllib.request.urlopen(http_address)
    out = os.path.join(genome_path, out_name + ".fa")
    tmp_out = os.path.join(genome_path, out_name + "tmp.fa")
    out_gz = out + ".gz"
    with open(out_gz,'wb') as outF:
        outF.write(opened.read())
    gf = gzip.open(out_gz)
    new_out = open(tmp_out,'wb')
    new_out.write(gf.read())
    gf.close()
    os.remove(out_gz)
    new_out.close()
    split_by_N(tmp_out, out, script)
    return out

"""
Given the created maps and the old config files, creates the required files and new config
"""
def write_config(otu_genome_map, no_samples, script, out_path_genomes):
    out_path = "./"
    genome_to_id = os.path.join(out_path, "genome_to_id.tsv")
    #config.set('community0','id_to_genome_file', genome_to_id)
    metadata = os.path.join(out_path, "metadata.tsv")
    with open(metadata,'w') as md:
        md.write("genome_ID\tOTU\tNCBI_ID\tnovelty_category\n") # write header
    #config.set('community0','metadata',metadata)
    #no_samples = int(config.get("Main","number_of_samples"))
    abundances = [os.path.join(out_path,"abundance_%s.tsv" % i) for i in range(no_samples)]
    #log.info("Downloading %s genomes" % len(otu_genome_map))
    
    #create_path = os.path.join(out_path_genomes,"genomes")
    create_path = out_path_genomes

    if not os.path.exists(create_path):
        os.makedirs(create_path)
    for otu in otu_genome_map:
        taxid, genome_id, path, novelty, curr_abundances = otu_genome_map[otu]
        counter = 0
        while counter < 10:
            try:
                if path.startswith('http') or path.startswith('ftp'):
                    genome_path = download_genome(path, out_path_genomes, script)
                else:
                    out_name = path.rstrip().split('/')[-1]
                    genome_path = os.path.join(create_path, out_name)
                    shutil.copy2(path, genome_path)
                break
            except Exception as e:
                counter += 1
                #log.error("Caught exception %s while moving/downloading genomes" % repr(e))
                if counter >= 10:
                    raise ValueError("Caught exception %s while moving/downloading genomes" % repr(e))
        if counter == 10:
#            log.error("Genome %s (from %s, path %s) could not be downloaded after 10 tries, check your connection settings" % (otu, genome_id, path))
            raise ValueError("Genome %s (from %s, path %s) could not be downloaded after 10 tries, check your connection settings" % (otu, genome_id, path))
        with open(genome_to_id,'a+') as gid:
            gid.write("%s\t%s\n" % (otu, genome_path))
        with open(metadata,'a+') as md:
            md.write("%s\t%s\t%s\t%s\n" % (otu,taxid,genome_id,novelty))
        i = 0
        for abundance in abundances:
            with open(abundance, 'a+') as ab:
                ab.write("%s\t%s\n" % (otu,curr_abundances[i]))
            i += 1


def str2bool(x):
    return str(x).lower() in ("1", "true", "yes", "y")


def set_ranks(ncbi):
    global RANKS

    # Get all possible taxonomic ranks
    conn = sqlite3.connect(ncbi.dbfile)
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT rank FROM species WHERE rank != ''")
    all_ranks = {row[0] for row in cursor.fetchall()}
    conn.close()

    if 'acellular root' in all_ranks:
        RANKS = NEW_RANKS
    else:
        RANKS = LEGACY_RANKS


def main(biom_profile, no_samples, reference_genomes, seed, mu, sigma,
         min_strains,  max_strains, debug, no_replace, fill_up, script,
         genomes_out_dir, additional_references, additional_genomes_quality_file, ncbi_taxdump_file,
         max_rank, prioritize_additional_genomes, additional_genomes_max_contamination,
         additional_genomes_min_completeness, additional_genomes_max_num_contigs):

    if additional_references == "None":
        additional_references = None
    if additional_genomes_quality_file == "None":
        additional_genomes_quality_file = None

    np_rand.seed(seed)

    ncbi = NCBITaxa(taxdump_file=ncbi_taxdump_file)
    set_ranks(ncbi)

    tax_profile = read_taxonomic_profile(biom_profile, no_samples)
    genomes_map, total_genomes = read_genomes_list(reference_genomes, additional_references,
                                                   additional_genomes_quality_file,
                                                   additional_genomes_max_contamination,
                                                   additional_genomes_min_completeness,
                                                   additional_genomes_max_num_contigs)
    per_rank_map = get_genomes_per_rank(genomes_map, ncbi)
    otu_genome_map, unmatched_otus = map_otus_to_genomes(tax_profile, per_rank_map, mu, sigma,
                                                         min_strains, max_strains, no_replace,
                                                         total_genomes, max_rank,
                                                         prioritize_additional_genomes, ncbi)

    if fill_up and len(unmatched_otus) > 0:
        otu_genome_map = fill_up_genomes(otu_genome_map, unmatched_otus, per_rank_map, tax_profile,
                                         debug, prioritize_additional_genomes)

    write_config(otu_genome_map, no_samples, script, genomes_out_dir)


if __name__ == "__main__":
    main(biom_profile = sys.argv[1],
        no_samples = int(sys.argv[2]),
        reference_genomes = sys.argv[3],
        seed = int(sys.argv[4]),
        mu = int(sys.argv[5]),
        sigma = int(sys.argv[6]),
        min_strains = int(sys.argv[7]),
        max_strains = int(sys.argv[8]),
        debug = str2bool(sys.argv[9]),
        no_replace = str2bool(sys.argv[10]),
        fill_up = str2bool(sys.argv[11]),
        script = sys.argv[12],
        genomes_out_dir = sys.argv[13],
        additional_references = sys.argv[14],
        additional_genomes_quality_file = sys.argv[15],
        ncbi_taxdump_file = sys.argv[16],
        max_rank = sys.argv[17],
        prioritize_additional_genomes = str2bool(sys.argv[18]),
        additional_genomes_max_contamination = float(sys.argv[19]),
        additional_genomes_min_completeness = float(sys.argv[20]),
        additional_genomes_max_num_contigs = int(float(sys.argv[21])))
