import sys
import os
import urllib
import gzip
import biom
import shutil
from numpy import random as np_rand
from ete3 import NCBITaxa
from scripts.loggingwrapper import LoggingWrapper as logger
from configparser import ConfigParser

def run_patch(): # patching ete3 version
    try:
        import ast
        import inspect
        import sys
        _log = logger()
        _log.info("Patching NCBITaxa's base methods. For reason, see https://github.com/etetoolkit/ete/issues/469.\n")
        code_to_patch = """db.execute("INSERT INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname))"""
        patched_code = """db.execute("INSERT OR REPLACE INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname))"""
        ncbiquery = sys.modules[NCBITaxa.__module__]
        lines_code = [x.replace(code_to_patch, patched_code)
                      for x in inspect.getsourcelines(ncbiquery.upload_data)[0]]
        # Insert info message to see if patch is really applied
        lines_code.insert(1, "    print('\\nIf this message shows, then the patch is successful!')\n")
        # Insert external import and constants since only this function is patched and recompiled
        lines_code.insert(1, "    import os, sqlite3, sys\n")
        lines_code.insert(1, "    DB_VERSION = 2\n")
        lines_code = "".join(lines_code)
        # Compile and apply the patch
        ast_tree = ast.parse(lines_code)
        patched_function = compile(ast_tree, "<string>", mode="exec")
        mod_dummy = {}
        exec(patched_function, mod_dummy)
        ncbiquery.upload_data = mod_dummy["upload_data"]
    except Exception as e:
        _log.info(e)
        _log.info("Patching failed, current taxonomy data downloaded from FTP may be failed to update with ETE3!")
    finally:
        _log.info("Patch finished.")
        _log = None

run_patch()
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()
RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
MAX_RANK = 'family'
_log = None

"""
Reads a BIOM file and creates map of OTU: lineage, abundance
BIOM file format needs to have a taxonomy field in metadata which contains the taxonomy in the format:
RANK__SCINAME; LOWERRANK_LOWERSCINAME
"""
def read_taxonomic_profile(biom_profile, config, no_samples = None):
    table = biom.load_table(biom_profile)
    ids = table.ids(axis="observation")
    samples = table.ids()

    if no_samples is None:
        no_samples = len(samples)

    if no_samples is not None and no_samples != len(samples) and no_samples != 1:
        _log.warning("Number of samples (%s) does not match number of samples in biom file (%s)" % (no_samples, len(samples)))
        if no_samples > len(samples):
            no_samples = len(samples)
        _log.warning("Using the first %s samples" % no_samples)

    config.set("Main", "number_of_samples", str(no_samples))
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

"""
Reads list of available genomes in the (tsv) format:
NCBI_ID Scientific_Name ftp_path
Additional files might be provided with:
NCBI_ID Scientific_Name genome_path novelty_category
were path might either be online or offline/local
"""
def read_genomes_list(genomes_path, additional_file = None):
    genomes_map = {}
    total_genomes = 0
    if additional_file is not None:
        with open(additional_file,'r') as add:
            for line in add:
                ncbi_id, sci_name, path, novelty = line.strip().split('\t')
                if ncbi_id in genomes_map:
                    genomes_map[ncbi_id][1].append(path)
                else:
                    genomes_map[ncbi_id] = (sci_name, [path], novelty) # this might not be a http path
                total_genomes += 1
    with open(genomes_path,'r') as genomes:
        for line in genomes:
            ncbi_id, sci_name, ftp = line.strip().split('\t')
            http = ftp.replace("ftp://","http://") # not using ftp address but http (proxies)
            if ncbi_id in genomes_map:
                genomes_map[ncbi_id][1].append(http)
            else:
                genomes_map[ncbi_id] = (sci_name, [http], 'known_strain') # sci_name is always the same for same taxid (?)
            total_genomes += 1
    return genomes_map, total_genomes

"""
Given all available genomes, creates a map sorted by ranks of available genomes on that particular rank, ordered by their ncbi ids
"""
def get_genomes_per_rank(genomes_map, ranks, max_rank):
    per_rank_map = {}
    for rank in ranks:
        per_rank_map[rank] = {}
    for genome in genomes_map:
        try:
            lineage = ncbi.get_lineage(genome) # this might contain some others ranks than ranks
            ranks_lin = ncbi.get_rank(lineage)
            for tax_id in lineage: # go over the lineage
                try:
                    check_rank = ranks_lin[tax_id]
                except KeyError:
                    continue
                if check_rank in per_rank_map: # if we are a legal rank
                    rank_map = per_rank_map[ranks_lin[tax_id]]
                    if tax_id in rank_map: # tax id already has a genome
                        for strain in genomes_map[genome][1]:
                            rank_map[tax_id].append((strain,genome)) # add http address
                    else:
                        rank_map[tax_id] = []
                        for strain in genomes_map[genome][1]:
                            rank_map[tax_id].append((strain,genome)) # add http address
        except ValueError as e:
           _log.warning(e)
    return per_rank_map

"""
Given a BIOM lineage, create a NCBI tax id lineage
"""
def transform_lineage(lineage, ranks, max_rank):
    new_lineage = []
    for member in lineage:
        name = member.split("__")[-1] # name is on the right hand side
        if len(name) == 0:
            continue
        mapping = ncbi.get_name_translator([name])
        if name in mapping:
            taxid = mapping[name][0]
            if ncbi.get_rank([taxid])[taxid] in ranks:
                new_lineage.append(taxid) # should contain only one element
        else:
            name = name.split()[0]
            if name in mapping:
                taxid = mapping[name][0]
                if ncbi.get_rank([taxid])[taxid] in ranks:
                    new_lineage.append(taxid) # retry if space in name destroys ID
    return new_lineage[::-1] # invert list, so lowest rank appears first (last in BIOM)

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

def all_genomes(per_rank_map):
    added_genomes = set()
    for rank in per_rank_map:
        for taxid in per_rank_map[rank]:
            for path, genome_id in per_rank_map[rank][taxid]:
                if path not in added_genomes:
                    added_genomes.add(path)
    _log.info(len(added_genomes))

"""
Given the OTU to lineage/abundances map and the genomes to lineage map, create map otu: taxid, genome, abundances
"""
def map_otus_to_genomes(profile, per_rank_map, ranks, max_rank, mu, sigma, max_strains, debug, no_replace, max_genomes):
    unmatched_otus = []
    otu_genome_map = {}
    warnings = []
    sorted_otus = sort_by_abundance(profile)
    genome_set_size = 0
    for otu in sorted_otus:
        if genome_set_size >= max_genomes and no_replace: #cancel if no genomes are available anymore
            break
        lin, abundances = profile[otu]
        lineage = transform_lineage(lin, ranks, max_rank)
        if len(lineage) == 0:
            warnings.append("No matching NCBI ID for otu %s, scientific name %s" % (otu, lin[-1].split("__")[-1]))
            unmatched_otus.append(otu)
        lineage_ranks = ncbi.get_rank(lineage)
        for tax_id in lineage: # lineage sorted ascending
            rank = lineage_ranks[tax_id]
            if ranks.index(rank) > ranks.index(max_rank):
                warnings.append("Rank %s of OTU %s too high, no matching genomes found" % (rank, otu))
                warnings.append("Full lineage was %s, mapped from BIOM lineage %s" % (lineage, lin))
                unmatched_otus.append(otu)
                break
            genomes = per_rank_map[rank]
            if tax_id not in genomes:
                warnings.append("For OTU %s no genomes have been found on rank %s with ID %s" % (otu, rank, tax_id))
                continue # warning will appear later if rank is too high
            available_genomes = genomes[tax_id]
            strains_to_draw = max((np_rand.geometric(2./max_strains) % max_strains),1)
            if len(available_genomes) >= strains_to_draw:
                used_indices = np_rand.choice(len(available_genomes),strains_to_draw,replace=False)
                used_genomes = set([available_genomes[i] for i in used_indices])
            else:
                used_genomes = set(available_genomes) # if not enough genomes: use all
            genome_set_size += len(used_genomes) # how many genomes are used
            log_normal_vals = np_rand.lognormal(mu,sigma, len(used_genomes))
            sum_log_normal = sum(log_normal_vals)
            i = 0
            for path, genome_id in used_genomes:
                otu_id = otu + "." + str(i)
                otu_genome_map[otu_id] = (tax_id, genome_id, path, []) # taxid, genomeid, http path, abundances per sample
                relative_abundance = log_normal_vals[i]/sum_log_normal
                i += 1
                for abundance in abundances: # calculate abundance per sample
                    current_abundance = relative_abundance * abundance
                    otu_genome_map[otu_id][-1].append(current_abundance)
                if (no_replace): # sampling without replacement:
                    for new_rank in per_rank_map:
                        for taxid in per_rank_map[new_rank]:
                            if (path, genome_id) in per_rank_map[new_rank][taxid]:
                                per_rank_map[new_rank][taxid].remove((path,genome_id))
            break # genome(s) found: we can break
    if len(warnings) > 0:
        _log.warning("Some OTUs could not be mapped")
        if debug:
            for warning in warnings:
                _log.warning(warning)
    return otu_genome_map, unmatched_otus, per_rank_map


"""
Take fasta input file and split by any N occurence (and remove Ns)
"""
def split_by_N(fasta_path, out_path):
    os.system("scripts/split_fasta.pl %s %s" % (fasta_path, out_path))
    os.remove(fasta_path)

"""
Downloads the given genome and returns the out path
"""
def download_genome(genome, out_path):
    genome_path = os.path.join(out_path,"genomes")
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
    split_by_N(tmp_out, out)
    return out

"""
Given the created maps and the old config files, creates the required files and new config
"""
def write_config(otu_genome_map, genomes_map, out_path, config):
    genome_to_id = os.path.join(out_path, "genome_to_id.tsv")
    config.set('community0','id_to_genome_file', genome_to_id)
    metadata = os.path.join(out_path, "metadata.tsv")
    with open(metadata,'w') as md:
        md.write("genome_ID\tOTU\tNCBI_ID\tnovelty_category\n") # write header
    config.set('community0','metadata',metadata)
    no_samples = int(config.get("Main","number_of_samples"))
    abundances = [os.path.join(out_path,"abundance%s.tsv" % i) for i in range(no_samples)]
    _log.info("Downloading %s genomes" % len(otu_genome_map))
    
    create_path = os.path.join(out_path,"genomes")
    if not os.path.exists(create_path):
        os.makedirs(create_path)
    for otu in otu_genome_map:
        taxid, genome_id, path, curr_abundances = otu_genome_map[otu]
        counter = 0
        while counter < 10:
            try:
                if path.startswith('http') or path.startswith('ftp'):
                    genome_path = download_genome(path, out_path)
                else:
                    out_name = path.rstrip().split('/')[-1]
                    genome_path = os.path.join(create_path, out_name)
                    shutil.copy2(path, genome_path)
                break
            except Exception as e:
                error = e
                counter += 1
                _log.error("Caught exception %s while moving/downloading genomes" % repr(e))
        if counter == 10:
            _log.error("Genome %s (from %s, path %s) could not be downloaded after 10 tries, check your connection settings" % (otu, genome_id, path))
        with open(genome_to_id,'a+') as gid:
            gid.write("%s\t%s\n" % (otu, genome_path))
        with open(metadata,'a+') as md:
            novelty = genomes_map[genome_id][-1]
            md.write("%s\t%s\t%s\t%s\n" % (otu,taxid,genome_id,novelty))
        i = 0
        for abundance in abundances:
            with open(abundance, 'a+') as ab:
                ab.write("%s\t%s\n" % (otu,curr_abundances[i]))
            i += 1
    abundance_files = ""
    for abundance in abundances[:-1]:
        abundance_files += abundance
        abundance_files += ","
    abundance_files += abundances[-1] # write csv of abundance files
    config.set("Main", 'distribution_file_paths', abundance_files)
    config.set("community0", "num_real_genomes", str(len(otu_genome_map)))
    config.set("community0", "genomes_total", str(len(otu_genome_map)))

    cfg_path = os.path.join(out_path, "config.ini")
    with open(cfg_path, 'w+') as cfg:
        config.write(cfg)
    return cfg_path

def fill_up_genomes(otu_genome_map, unmatched_otus, per_rank_map, tax_profile, debug):
    genomes = {}
    added_genomes = set()
    for rank in per_rank_map:
        for taxid in per_rank_map[rank]:
            genomes[taxid] = []
            for path, genome_id in per_rank_map[rank][taxid]:
                if path not in added_genomes:
                    genomes[taxid].append((path, genome_id))
                    added_genomes.add(path)
    otu_indices = np_rand.choice(len(unmatched_otus),len(unmatched_otus),replace=False)
    i = 0
    set_all = False
    for tax_id in genomes:
        for path, genome_id in genomes[tax_id]:
            curr_otu = unmatched_otus[otu_indices[i]] #so we choose a random genome
            lineage, abundances = tax_profile[curr_otu]
            lin = transform_lineage(lineage, RANKS, MAX_RANK)
            otu_genome_map[curr_otu] = (tax_id, genome_id, path, abundances)
            if debug:
                _log.warning("Filling up OTU %s (mapped tax id: %s) to genome with tax id %s" % (curr_otu, lin[0], tax_id))
            i += 1
            if (i >= len(unmatched_otus) or i >= len(added_genomes)):
                set_all = True
                break
        if (set_all):
            break
    return otu_genome_map

def generate_input(args):
    global _log
    _log = logger(verbose = args.debug)
    np_rand.seed(args.seed)
    #MAX_RANK = args.maxrank
    config = ConfigParser()
    config.read(args.config)
    try:
        max_strains = int(config.get("Main", max_strains_per_otu))
    except:
        max_strains = 3 # no max_strains have been set for this community - use cami value
        _log.warning("Max strains per OTU not set, using default (3)")
    try:
        mu = int(config.get("Main", "log_mu"))
        sigma = int(config.get("Main", "log_sigma"))
    except:
        mu = 1
        sigma = 2 # this aint particularily beatiful
        _log.warning("Mu and sigma have not been set, using defaults (1,2)") #TODO 
    tax_profile = read_taxonomic_profile(args.profile, config, args.samples)
    genomes_map, total_genomes = read_genomes_list(args.reference_genomes, args.additional_references)
    per_rank_map = get_genomes_per_rank(genomes_map, RANKS, MAX_RANK)
    otu_genome_map, unmatched_otus, per_rank_map = map_otus_to_genomes(tax_profile, per_rank_map, RANKS, MAX_RANK, mu, sigma, max_strains, args.debug, args.no_replace, total_genomes)
    if (args.fill_up and len(unmatched_otus) > 0):
        otu_genome_map = fill_up_genomes(otu_genome_map, unmatched_otus, per_rank_map, tax_profile, args.debug)
    cfg_path = write_config(otu_genome_map, genomes_map, args.o, config)
    _log.info("Community design finished")
    _log = None
    return cfg_path

