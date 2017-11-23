import sys
import os
from ftplib import FTP 
import gzip
import random
import biom #TODO required software
from numpy import random as np_rand
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy
from scripts.Validator.validator import Validator
from scripts.loggingwrapper import LoggingWrapper as logger
try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser  # ver. < 3.0

"""
Given a 16S-profile (currently only in CAMI format), downloads all closest relative genomes and creates abundances
"""

# strain inclusion?
RANKS=['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
#map BIOM ranks to CAMI ranks
BIOM_RANKS={'s':0,'g':1,'f':2,'o':3,'c':4,'p':5,'k':6}
THRESHOLD="family" #level up to which related genomes are to be found
_log = None

"""
Reads a biom (from e.g. QIIME) profile and transforms it so it can be used for the pipeline
"""
def transform_profile(biom_profile, no_samples, taxonomy):
    try:
        table = biom.load_table(biom_profile)
    except:
        try:
            return read_profile(biom_profile) # file is not a biom file: CAMI format
        except:
            _log.error("Incorrect file format of input profile")
            return
    ids = table.ids(axis="observation")
    samples = table.ids() # the samples' ids of the biom file
        
    profile = {} 
    warnings_rank = [] # collect all warnings
    warnings_sciname = [] # if scientific name wasnt found
    if no_samples is not None and no_samples != len(samples) and no_samples != 1: # no. samples not equal to samples in biom file, simulate using only the first sample
        _log.warning("Number of samples in biom file does not match number of samples in biom file, using first biom sample for simulation")
        no_samples = 1
    elif no_samples is None:
        no_samples = len(samples)
    for id in ids:
        lineage = table.metadata(id,axis="observation") # retrieving lineage
        if lineage is None: 
            lineage = id.split(";") # in prepared biom files the id is already the taxonomy TODO (might need to split by | or other char)
        else:
            lineage = lineage["taxonomy"] # "original" biom file stores taxonomy in metadata/taxonomy
            try:
                lineage = lineage.split(";") # sometimes this is still needed, grr
            except:
                pass
        abundances = []
        for sample in samples[:no_samples]:
            metadata = []
            abundances.append(table.get_value_by_ids(id,sample))
        ncbi_id, tax_path, sci_name = map_to_ncbi_id(lineage, taxonomy)
        if ncbi_id is None:
            if not tax_path: # rank was too high
                warnings_rank.append((RANKS[BIOM_RANKS[sci_name[0]]],sci_name[1]))
            else:
                warnings_sciname.append(sci_name)
            continue # do not add empty hits
        profile[id] = (ncbi_id, tax_path, abundances)
    
    if len(warnings_rank):
        _log.warning("Some genomes had a too high rank and were omitted")
        for warning in warnings_rank:
            _log.info("Rank (%s) of genome %s was too high" % warning)
    if len(warnings_sciname):
        _log.warning("Some scientific names were not found and omitted")
        for warning in warnings_sciname:
            _log.info("Scientific name %s did not match any in NCBI" % warning)
    return profile, no_samples

"""Given the biom lineage, calculate the NCBI lineage"""
def map_to_ncbi_id(lin, taxonomy):
    lineage = []
    for rank in lin: # only the lineage
        taxon = rank.strip().split("__") # split biom-string
        if len(taxon) != 2: # there is no name
            break
        if taxon[1] == '': 
            if BIOM_RANKS[taxon[0]] >= RANKS.index(THRESHOLD): # the rank is higher than the desired threshold
                return None, False, lineage[-1]
        lineage.append(taxon)
    sci_name = retrieve_scientific_name(lineage, False)
    sci_name.encode('ascii','ignore') # and hope that this does not break something
    sci_name = str(sci_name) # since it has been encoded this cast shouldnt fail
    ncbi_ids = taxonomy.get_taxids_by_scientific_name(sci_name, True)
    #ncbi_ids = taxonomy.get_taxids_by_scientific_name_wildcard(sci_name)
    if ncbi_ids is None:
        sci_name = retrieve_scientific_name(lineage, True)
        sci_name.encode('ascii','ignore') # and hope that this does not break something
        sci_name = str(sci_name) # since it has been encoded this cast shouldnt fail
        ncbi_ids = taxonomy.get_taxids_by_scientific_name(sci_name,True)
        #ncbi_ids = taxonomy.get_taxids_by_scientific_name_wildcard(sci_name)
        if ncbi_ids is None:
            return None, True, sci_name
    ncbi_id = ncbi_ids.pop() # TODO do not take first if more than one?
    tax_path = taxonomy.get_lineage_of_legal_ranks(ncbi_id)
    return ncbi_id, tax_path, sci_name

"""
Given the biom-lineage information retrieves the scientific name, which is composed of genus + species for species level, or just the information of the lowest rank otherwise
"""
def retrieve_scientific_name(lineage, try_other_format):
    name = ''
    for rank in lineage:
        sci_name = rank[1]
        if sci_name != '':
            name = sci_name # lowest set rank is used
        if try_other_format and rank[0] == 's' and rank[1] != '': 
            name += " " + rank[1] # species name starts with genus name (e.g. genus Escherichia, species coli)
    return name

"""
original code in the profiling-evaluation-biobox, reads a file in cami profiling format and extracts relevant information (taxids/tax path/relative abundance/genome rank in the taxonomy)
adpated from the profile evaluation biobox, extendeded by the following: We only check for the species tax ids, original genomes for the higher ranks will be checked later on
"""
def read_profile(file_path):
    if not isinstance(file_path, basestring):
        _log.error("file_path is invalid: %s" % file_path)
    if isinstance(file_path, str) and not os.path.isfile(file_path):
        _log.error("16S profile not found in: %s" % file_path)
        raise Exception("File not found")
    # check whether profile is biom or cami format
    
    profile = {}
    with open(file_path, 'r') as read_handler:
        for line in read_handler:
            line = line.rstrip()
            if len(line) == 0:
                continue  # skip blank lines
            if line.lower().split(':')[0] == '@ranks':
                ranks = line.strip().split(':')[1].split('|')
                continue
            if line[0] in ['@', '#']:
                continue  # skip comment or header
            profile_info  = line.split('\t')
            taxid = profile_info[0]
            rank = profile_info[1]
            taxpath = profile_info[2]
            taxpath_sciname = profile_info[3]
            sciname = taxpath_sciname.split('|')[-1] # deepest scientific name
            weight = float(profile_info[4])
            if weight == 0:
                # Ignore zero weighted taxIDs
                continue
            if (rank == 'species'): # only search for genomes on species level, make variable
                profile[sciname] = (taxid, taxpath, [weight]) # cami format currently only supports single sample
    return profile 

"""
given the list of full genomes available from NCBI, create a mapping with the relevant data (ncbi id/scientific name/ftp address of full genomes)
file_path is the path to a file created from the NCBI assembly summary: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
This file is filtered such that only assemblies specified with "Complete genome" are left and the only three columns are
NCBI tax id \t scientific name \t ftp address
scientific name is for debugging, the ftp address the address of the corresponding full genome for downloading
"""
def read_genome_list(file_path):
    assert isinstance(file_path, basestring)
    if isinstance(file_path, str) and not os.path.isfile(file_path):
        _log.error("Reference genome list not found in: %s" % file_path)
        raise Exception("File not found")
    tax_ids = list()
    sci_name = {}
    ftp_address = {}
    with open(file_path,'r') as full_genomes:
        for line in full_genomes:
            if len(line) == 0:
                continue
            temp = line.split('\t')
            if len(temp) < 3: # if the line split by tabs does not contain the three elements something is wrong
                continue
            tax_ids.append(temp[0])
            sci_name[temp[1]] = temp[0]
            ftp_address[temp[0]] = temp[2]
    return tax_ids, sci_name, ftp_address
"""
extends the list of tax ids to a list of list to include higher level taxonomic ranks
in the extended genome list the first level is just the list of all tax ids for which a complete genome is present
for all of the ranks starting from genus, ..., superkingdom a dictionary is created, mapping from a tax id on the
level of the current rank to the species' tax ids belonging to that higher rank tax id, i.e.
Assume the tax id of a given species is "83333" (Escherichia coli K-12), then "83333" is present on the first level, since 
a full genome is known.
Escherichia coli K-12 belong to the genus Escherichia with tax id "561", so there is a mapping:
"561" : ["83333"]
Additionally, assume there is a full genome of Escherichia albertii available (tax id "208962"), then the mapping would be
"561" : ["83333", "208962"]
"""
def extend_genome_list(tax_ids, tax):
    ids_per_rank = [tax_ids] # on the first rank no mapping
    for rank in RANKS[1:]:
        ids_per_rank.append(dict()) # create empty dicts
    for tax_id in tax_ids:
        try:
            lineage = tax.get_lineage_of_legal_ranks(tax_id, ranks = RANKS)
        except ValueError:
            continue # this means, the taxid was not found in the reference (reference mismatch)
        i = 1
        for rank_tax_id in lineage[1:]: # for the path up to the source/superkingdom
            if rank_tax_id in ids_per_rank[i]:# add the species id to the map, if the higher rank tax id is present
                ids_per_rank[i][rank_tax_id].append(tax_id) 
            else:
                ids_per_rank[i].update({rank_tax_id : [tax_id]})
            if rank == THRESHOLD: # only search up to this taxonomic rank
                break 
            i = i + 1 # continue with the next rank
    return ids_per_rank

""" 
Given the list of available full genomes and the species extracted from the profile, finds a mapping of profile genomes to available full genomes sequences
First, extends the list of reference tax ids, so the least common ancestors in the phylogeny can be found.
For all species tax ids in the profile, it is checked whether that tax id is a tax id of a known full genome.
If not, than the next rank is checked, i.e. genus. If there is a full genome with the same genus tax id like our genus tax id,
than one of these genomes is chosen as the "closest related" genome. 
"""
def map_to_full_genomes(ref_tax_ids, profile, tax, seed):
    extended_genome_list = extend_genome_list(ref_tax_ids,tax)
    to_download = dict() # ncbi id of genomes to download
    _log.info("Downloading genomes from NCBI")
    warnings = [] # warnings if no complete genomes are found
    for ids in profile:
        taxid = profile[ids][0]
        found_genome = False
        if taxid in extended_genome_list[0]: # a full genome with exact ncbi id is present
            to_download.update({ids : ([taxid], taxid)}) 
            found_genome = True
        else: # the exact genome is not present, go up the ranks
            try:
                lineage = tax.get_lineage_of_legal_ranks(taxid,ranks = RANKS)
            except ValueError: #tax ID was not found in reference data base
                _log.warning("Genome %s not found in reference, maybe your reference is deprecated?" % taxid)
                continue
            i = 0
            for higher_taxid in lineage: # rank is a number corresponding to the ranks defined in RANKS with species being the lowers (0)
                if higher_taxid is not None and higher_taxid in extended_genome_list[i]:
                    species_id = extended_genome_list[i][higher_taxid]
                    to_download.update({ids : (species_id, higher_taxid)}) #all matching genomes
                    found_genome = True
                    break #TODO add rank for debugging purposes (RANKS[i])
                if RANKS[i] == THRESHOLD:
                    break # we dont need to search on higher levels
                i += 1 # go to the next rank
        if not found_genome: # No reference genome up until THRESHOLD
            warnings.append(taxid)
    if len(warnings):
        _log.warning("Some NCBI IDs did not map to complete genomes")
        for warning in warnings:
            _log.info("No genome corresponding to ID %s found, omitted." % warning)
    return to_download

"""
Given the ftp address of a list of full genomes, download the corresponding genomes from ncbi
The ftp server address of ncbi is ftp.ncbi.nlm.nih.gov (make sure this did not change)
iterates over the list of genomes in the profile, retrieves the mapped full genome and its ftp address
based from this downloads the file.
The path contains more files, the sequence ends with _genomic.fna.gz
We might also download the _genomic.gff.gz for genes/evolution
"""
def download_genomes(list_of_genomes, ftp_list, out_path):
    metadata = dict() # create the metadata table (pathes to genomes)
    warnings = []
    ftp = FTP('ftp.ncbi.nlm.nih.gov') #reduce timeout?
    ftp.login() # anonymous login
    sample_path = os.path.join(out_path,"genomes")
    _log.info("Downloading %s genomes" % len(list_of_genomes))
    if not os.path.exists(sample_path):
        os.makedirs(sample_path)
    for genome_id in list_of_genomes:
        gen = list_of_genomes[genome_id][0]
        otu = list_of_genomes[genome_id][1]
        path = ftp_list[gen]
        split_path = path.split('/')
        cwd = "/" + "/".join(split_path[3:]).rstrip() # get /address/to/genome
        gen_name = split_path[-1].rstrip() # genome name is last in address
        to_dl = gen_name + "_genomic.fna.gz"
        out_name = os.path.join(sample_path,gen) + ".fa"  # out name is the ncbi id of the downloaded genome
        metadata.update({genome_id:(out_name, otu)})
        if (os.path.isfile(out_name)): # we already downloaded this genome
            continue
        out_name_gz = out_name + ".gz"
        try:
            ftp.cwd(cwd)
        except: # huh, lets try again
            ftp = FTP('ftp.ncbi.nlm.nih.gov') #reduce timeout?
            ftp.login() # anonymous login
            ftp.cwd(cwd)
        counter = 0
        while (counter < 10):
            try: 
                ftp.retrbinary("RETR %s" % to_dl, open(out_name_gz,'wb').write) #download genomes
                break
            except:
                counter += 1
        if (counter == 10):
            warnings.append("File %s could not be downloaded (Genome ID %s/NCBI ID %s" % (to_dl,genome_id,out_name))
            metadata[genome_id] = None
            continue
        gf = gzip.open(out_name_gz) 
        outF = open(out_name,'wb')
        outF.write(gf.read())
        gf.close()
        os.remove(out_name_gz) # remove the now unzipped archives
        outF.close()
    if len(warnings):
        log.warning("Downloading %s genomes failed, try running with --debug if this happends regularily" % len(warnings))
        for warning in warnings:
            log.debug(warning)
    return metadata

"""
Given the list of genomes and the profile, create an abundance table for the downloaded genomes
"""
def create_abundance_table(list_of_genomes, seed, config, profile):
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
    np_rand.seed(seed)
    if max_strains >= 2:
        strains_to_draw = max((np_rand.geometric(2./max_strains) % max_strains),1) # make sure we draw at least one
    else:
        strains_to_draw = 1
    # mean will be ~max_strains/2 and a minimum of 1 strain is drawn
    abundances = {}
    genome_selection = {}
    for elem in list_of_genomes:
        total_abundances = profile[elem][2] # profile has a taxid - weight map at pos 2
        mapped_genomes = list_of_genomes[elem][0]
        otu = list_of_genomes[elem][1]
        if len(mapped_genomes) == 0:
            _log.warning("All mapping genomes for OTU %s have been used" % otu)
            continue
        if len(mapped_genomes) >= strains_to_draw: # if more genomes mapped than needed, do a selection
            mapped_genomes = [mapped_genomes[x] for x in np_rand.choice(len(mapped_genomes), strains_to_draw)] # sample genomes    
        log_normal_vals = np_rand.lognormal(mu,sigma,len(mapped_genomes))
        sum_log_normal = sum(log_normal_vals)
        i = 0
        for g in mapped_genomes:
            relative_abundance = log_normal_vals[i]/sum_log_normal
            genome_id = elem + "." + str(i)
            i += 1
            for abundance in total_abundances:
                current_abundance = relative_abundance * abundance
                genome_selection[genome_id] = (g, otu)
                if genome_id not in abundances:
                    abundances[genome_id] = [current_abundance] # relative abundance in lognormal (%) times total_abundance (in profile)
                else:
                    abundances[genome_id].append(current_abundance)
        for used_genome in mapped_genomes: # sample without replacement (remove used)
            list_of_genomes[elem][0].remove(used_genome) 
    return abundances, genome_selection

# read args
def read_args(args):    
    genome_list = args.reference_genomes
    profile = args.profile
    tax_path = args.ncbi
    #download = args.dont_download_genomes
    seed = args.seed
    no_samples = args.samples
    out_path = os.path.join(args.o,'') #so we are sure it is a directory
    config = ConfigParser()
    config.read(args.config)
    tax = NcbiTaxonomy(tax_path)
    #return genome_list, profile, tax_path, download, seed, no_samples, out_path, config, tax
    return genome_list, profile, tax_path, seed, no_samples, out_path, config, tax

"""given all samples' profiles, downloads corresponding genomes and creates required tables
abundance: mapping from downloaded_genomes to their abundance
mapping: mapping from genome_id to list of all downloaded genomes and their otu
downloaded: mapping from genome_id to file path"""
def create_full_profiles(profile, tid, ftp, tax, seed, config, out_path):
    to_download = map_to_full_genomes(tid,profile,tax,seed) #returns id:[mapped_genomes],otu
    abundances, genome_selection = create_abundance_table(to_download,seed,config,profile)
    genomes = download_genomes(genome_selection,ftp,out_path)
    return genomes, abundances

"""creates and writes the files required for configuration"""
def create_configs(out_path, config, abundances, downloaded, nr_samples):
    filenames = []
    for i in xrange(nr_samples): 
        filename = os.path.join(out_path,"abundance%s.tsv" % i)
        with open(filename,'wb') as abundance_i:
            for genome in abundances:
                abundance = abundances[genome][i]
                if genome not in downloaded:
                    _log.warning("Genome with abundance %s was not downloaded" % abundance)
                    continue # this has not been downloaded
                abundance_i.write("%s\t%s\n" % (genome,abundance))
        filenames.append(filename)

    filename = os.path.join(out_path,"genome_to_id.tsv")
    with open(filename,'wb') as gpath:
        for genome_id in downloaded:
            tax_id = downloaded[genome_id][0]
            gpath.write("%s\t%s\n" % (genome_id,tax_id))
    config.set('community0','id_to_genome_file',filename)
   
    filename = os.path.join(out_path,"metadata.tsv")
    with open(filename,'wb') as metadata:
        metadata.write("genome_ID\tOTU\tNCBI_ID\tnovelty_category\n") # header
        for genome_id in downloaded:
            path_to_genome = downloaded[genome_id][0]
            ncbi_id = path_to_genome.rsplit("/",1)[-1].rsplit(".",1)[0] # split at path and then strip file ending
            otu = downloaded[genome_id][1]
            metadata.write("%s\t%s\t%s\t%s\n" % (genome_id,otu,ncbi_id,"new_strain")) #check multiple matchings
    config.set('community0','metadata',filename)
    
    config.set('community0','num_real_genomes',str(len(downloaded)))
    config.set('community0','genomes_total',str(len(downloaded)))
     # TODO what if strains should be simulated?
     # TODO error if genome_total and num_real are set but too small

    filename_list = "" # write all filenames as comma-separated list
    for filename in filenames[:-1]:
        filename_list += filename
        filename_list += ","
    filename_list += filenames[-1]
    
    config.set('Main', 'distribution_file_paths', filename_list) # distributions of all samples
    config.set('Main', 'number_of_samples', nr_samples)

    cfg_path = out_path + "config.ini"
    with open(cfg_path,'wb') as cfg:
        config.write(cfg)
    return cfg_path

"""
Given the reference genomes' sequences (path), an 16S profile, the path to the NCBI taxonomy and the output path,
downloads mapped genomes, creates an abundance table and all the further inputs which are needed downstream by the main pipeline.
If download is set to false, no genomes are downloaded and instead the reference genomes are expected to be in the out directory
The file name should then be out_path/taxID.fa.gz so it can be found
"""
#list of full genomes, input profile (CAMI format), taxonomy path, out directory
#returns number of genomes
def generate_input(args):
    #genome_list, profile, tax_path, download, seed, no_samples, out_path, config, tax = read_args(args)
    global _log
    _log = logger(verbose = args.debug)
    genome_list, profile, tax_path, seed, no_samples, out_path, config, tax = read_args(args)

    tax_ids, sci_names, ftp = read_genome_list(genome_list)
    
    profile, nr_samples = transform_profile(profile,args.samples,tax) # might be multiple ones if biom file
    
    downloaded, abundances = create_full_profiles(profile, tax_ids, ftp, tax, seed, config, out_path)

    res = create_configs(out_path, config, abundances, downloaded, nr_samples)
    _log = None
    return res

