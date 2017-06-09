import sys
import os
from ftplib import FTP 
import gzip
import random
import biom #TODO required software
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
log = logger(verbose=False)

"""
Reads a biom (from e.g. QIIME) profile and transforms it so it can be used for the pipeline
"""
def transform_profile(biom_profile, epsilon, no_samples, taxonomy):
	try:
		table = biom.load_table(biom_profile)
	except:
		return read_profile(biom_profile, epsilon) # file is not a biom file: CAMI format
	ids = table.ids(axis="observation")
	profiles = []
	samples = table.ids() # the samples' ids of the biom file
	i = 0
	for sample in samples:
		metadata = []
		log.info("Processing sample %s" % i)
		if no_samples is not None and no_samples != len(samples) and no_samples != 1 and i > 0: # no. samples not equal to samples in biom file, simulate using only the first sample
			log.warning("Number of samples in biom file does not match number of samples in biom file, using first biom sample for simulation")
			break
		#if i > no_samples: # simulate the first i samples from biom file if i < no_samples?
		#	break
		profile = ([],[],{},RANKS) # tax ids / tax paths / tax_id : weight map / ranks
		for id in ids:
			abundance = table.get_value_by_ids(id,sample)
			lineage = table.metadata(id,axis="observation") # retrieving lineage
			if lineage is None: 
				lineage = id.split(";") # in prepared biom files the id is already the taxonomy TODO
			else:
				lineage = lineage["taxonomy"] # "original" biom file stores taxonomy in metadata/taxonomy
			if abundance > 0: #only present strains should appear in profile
				if len(lineage) > 1: # if length of lineage is one then the strain cannot be assigned
					metadata.append((lineage,abundance))
		for elem in metadata:
			lineage = []
			for rank in elem[0]: # only the lineage
				tax = rank.split("__") # split biom-string
				if len(tax) != 2: # there is no name
					break
				if tax[1] == '': 
					if BIOM_RANKS[tax[0]] >= RANKS.index(THRESHOLD): # the rank is higher than the desired threshold
						log.warning("Rank (%s) of %s is too high, omitted." % (RANKS[BIOM_RANKS[lineage[-1][0]]],lineage[-1][1]))
						lineage = [] # skip this genome
					break #so we get the lowest set rank (assuming no rank is bypassed)
				lineage.append(tax)
			if lineage == []: # rank is too high, ignore
				continue
			sci_name = retrieve_scientific_name(lineage)
			sci_name.encode('ascii','ignore') # and hope that this does not break something
			sci_name = str(sci_name) # since it has been encoded this cast shouldnt fail
			ncbi_ids = taxonomy.get_taxids_by_scientific_name_wildcard(sci_name)# which one if there is more than one?
			# does this speed up/find less when not wildcard is used? TODO
			if ncbi_ids is None:
				log.warning("Scientific name %s does not correspond to NCBI id, omitted" % sci_name)
				continue
			ncbi_id = ncbi_ids.pop()#TODO
			tax_path = taxonomy.get_lineage_of_legal_ranks(ncbi_id) 
			weight = elem[1]
			profile[0].append(ncbi_id)
			profile[1].append(tax_path)
			if ncbi_id in profile[2]:
				profile[2][ncbi_id] += weight
			else:
				profile[2][ncbi_id] = weight
		profiles.append(profile)
	return profiles

"""
Given the biom-lineage information retrieves the scientific name, which is composed of genus + species for species level, or just the information of the lowest rank otherwise
"""
def retrieve_scientific_name(lineage):
	name = ''
	for rank in lineage:
		sci_name = rank[1]
		if sci_name != '':
			name = sci_name # lowest set rank is used
		#if rank[0] == 's' and rank[1] != '': #TODO check how this is working
		#	name += " " + rank[1] # species name starts with genus name (e.g. genus Escherichia, species coli)
	return name

"""
original code in the profiling-evaluation-biobox, reads a file in cami profiling format and extracts relevant information (taxids/tax path/relative abundance/genome rank in the taxonomy)
adpated from the profile evaluation biobox, extendeded by the following: We only check for the species tax ids, original genomes for the higher ranks will be checked later on
"""
def read_profile(file_path, epsilon):
	assert isinstance(file_path, basestring)
	if isinstance(file_path, str) and not os.path.isfile(file_path):
		log.error("16S profile not found in: %s" % file_path)
		raise Exception("File not found")
	assert epsilon is None or isinstance(epsilon, (float, int, long))
	# check whether profile is biom or cami format
	
	tax_path = list()
	tax_ids = list()
	weights = dict()
	ranks = []
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
			temp_split = line.split('\t')
			weight = float(temp_split[4])
			if epsilon is not None and weight < epsilon:
				# Ignore taxIDs below cutoff
				continue
			if weight == 0:
				# Ignore zero weighted taxIDs
				continue
			if (temp_split[1] == 'species'): # only search for genomes on species level
				tax_path.append(temp_split[2])  # add the whole taxpath
				tax_ids.append(temp_split[0])  # just terminal tax ID
				weights[temp_split[0]] = weight  # the associated weight
	return [(tax_ids, tax_path, weights, ranks)] # since in biom this might be multiple ones

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
		log.error("Reference genome list not found in: %s" % file_path)
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
def map_to_full_genomes(ref_tax_ids, profile_tax_ids, tax, sample, seed):
	random.seed(seed)
	gen_map = {}
	extended_genome_list = extend_genome_list(ref_tax_ids,tax)
	to_download = dict() # ncbi id of genomes to download
	log.info("Downloading genomes from NCBI")
	for taxid in profile_tax_ids[0]:
		found_genome = False
		if taxid in extended_genome_list[0]: # a full genome with exact ncbi id is present
			to_download.update(taxid = taxid) # if profile contains strains this might cause overwrites TODO
			found_genome = True
		else: # the exact genome is not present, go up the ranks
			try:
				lineage = tax.get_lineage_of_legal_ranks(taxid,ranks = RANKS)
			except ValueError: #tax ID was not found in reference data base
				log.warning("Genome %s not found in reference, maybe your reference is deprecated?" % taxid)
				continue
			i = 0
			for higher_taxid in lineage: # rank is a number corresponding to the ranks defined in RANKS with species being the lowers (0)
				if higher_taxid in extended_genome_list[i]:
					species_id = extended_genome_list[i][higher_taxid]
					to_download.update({taxid : species_id[random.randint(0,len(species_id) - 1)]}) #randomly select one of the mapped genomes TODO 
					found_genome = True
					break #TODO add rank for debugging purposes (RANKS[i])
				if RANKS[i] == THRESHOLD:
					break # we dont need to search on higher levels
				i = i + 1 # go to the next rank
		if not found_genome: # No reference genome up until THRESHOLD
			log.warning("No genome corresponding to ID %s found, omitted." % taxid)
	return to_download

"""
Given the ftp address of a list of full genomes, download the corresponding genomes from ncbi
The ftp server address of ncbi is ftp.ncbi.nlm.nih.gov (make sure this did not change)
iterates over the list of genomes in the profile, retrieves the mapped full genome and its ftp address
based from this downloads the file.
The path contains more files, the sequence ends with _genomic.fna.gz
We might also download the _genomic.gff.gz for genes/evolution
Also note that, if by chance multiple original genomes mapped to the same reference genome, this will get downloaded multiple times,
but should only appear once in the out directory.
"""
def download_genomes(list_of_genomes, ftp_list, out_path, sample):
	metadata = dict() # create the metadata table (pathes to genomes)
	ftp = FTP('ftp.ncbi.nlm.nih.gov') 
	ftp.login() # anonymous login
	sample_path = os.path.join(out_path,"sample%s" % sample) # extra folder for every sample
	log.info("Downloading %s genomes" % len(list_of_genomes))
	if not os.path.exists(sample_path):
		os.makedirs(sample_path)
	out_path = os.path.join(sample_path,"genomes") # extra folder for downloaded genomes
	if not os.path.exists(out_path):
		os.makedirs(out_path)
	for elem in list_of_genomes:
		path = ftp_list[list_of_genomes[elem]]
		split_path = path.split('/')
		cwd = "/" + "/".join(split_path[3:]).rstrip() # get /address/to/genome
		gen_name = split_path[-1].rstrip() # genome name is last in address
		to_dl = gen_name + "_genomic.fna.gz"
		out_name = os.path.join(out_path,list_of_genomes[elem]) + ".fa"  # out name is the ncbi id of the downloaded genome
		out_name_gz = out_name + ".gz"
		metadata.update({list_of_genomes[elem]:out_name})
		ftp.cwd(cwd)
		ftp.retrbinary("RETR %s" % to_dl, open(out_name_gz,'wb').write) #download genomes
		gf = gzip.open(out_name_gz) 
		outF = open(out_name,'wb')
		outF.write(gf.read())
		gf.close()
		os.remove(out_name_gz) # remove the now unzipped archives
		outF.close()
	return metadata

"""
Given the list of genomes and the profile, create an abundance table for the downloaded genomes
"""
#TODO log-normal distribution within a genome
def create_abundance_table(list_of_genomes, profile):
	abundance = {}
	for elem in list_of_genomes:
		ab = float(profile[2][elem]) # profile has a taxid - weight map at pos 2
		ab = int(ab*1000) # just so we get nicer numbers
		abundance.update({list_of_genomes[elem]:ab})
	return abundance

# read args
def read_args(args):	
	genome_list = args.reference_genomes
	profile = args.profile
	tax_path = args.ncbi
	download = args.dont_download_genomes
	seed = args.seed
	no_samples = args.samples
	out_path = os.path.join(args.o,'') #so we are sure it is a directory
	config = ConfigParser()
	config.read(args.config)
	tax = NcbiTaxonomy(tax_path)
	return genome_list, profile, tax_path, download, seed, no_samples, out_path, config, tax

"""given all samples' profiles, downloads corresponding genomes and creates required tables"""
def create_full_profiles(profiles, tid, ftp, tax, seed, download, out_path):
	i = 0
	abundances = []
	downloaded = []
	mapping = []
	for profile in profiles:
		mapping.append(map_to_full_genomes(tid,profile,tax,i,seed))
		to_dl = mapping[i]

		if not download:
			downloaded.append({})
			for gen in to_dl:
				genome = to_dl[gen]
				downloaded[i].update({genome:os.path.join(out_path,"sample%s" % i,"genomes", to_dl[gen]) + ".fa" % i})
		else:
			downloaded.append(download_genomes(to_dl,ftp,out_path,i))
		abundances.append(create_abundance_table(to_dl,profile))
		i += 1
	return downloaded, abundances, mapping, i

"""creates and writes the files required for configuration"""
def create_configs(i, out_path, config, abundances, downloaded, mapping):
	numg = 0
	for k in xrange(i): # number of samples TODO samples vs #profiles!
		sample_path = out_path + "sample%s/" % k
		current_community = 'community%s' % k
		if current_community not in config.sections():
			config.add_section(current_community)
			section_items = config.items('community0')
			for item in section_items:
				config.set(current_community,item[0],item[1]) # set the other values for the communities like in the first community
		#TODO this can also be customized

		with open(sample_path + "abundance.tsv",'wb') as abundance:
			for gen in mapping[k]:
				genome = mapping[k][gen]
				if abundances[k][genome] == 0: # abundance is too low, do not simulate reads
					continue
				abundance.write("%s\t%s\n" % (gen,abundances[k][genome]))
		filename = sample_path + "abundance.tsv"
		config.set(current_community, 'distribution_file_paths',filename)

		with open(sample_path + "genome_to_id.tsv",'wb') as gpath:
			for gen in mapping[k]:
				genome = mapping[k][gen]
				gpath.write("%s\t%s\n" % (gen,downloaded[k][genome]))
		filename = sample_path + "genome_to_id.tsv"
		config.set(current_community,'id_to_genome_file',filename)
		
		with open(sample_path + "metadata.tsv",'wb') as metadata:
			metadata.write("genome_ID\tOTU\tNCBI_ID\tnovelty_category\n") # header
			otu = 0 #for OTU assignment (every species gets its own OTU here) TODO
			for gen in mapping[k]:
				metadata.write("%s\t%s\t%s\t%s\n" % (gen,otu,mapping[k][gen],"new_strain"))
				otu = otu + 1
		filename = sample_path + "metadata.tsv"
		config.set(current_community,'metadata',filename)
		
		config.set(current_community,'genomes_total',str(len(downloaded[k]))) # TODO what if strains should be simulated
		
		numg += len(downloaded[k])
	cfg_path = out_path + "config.ini"
	with open(cfg_path,'wb') as cfg:
		config.write(cfg)
	return numg, cfg_path

"""
Given the reference genomes' sequences (path), an 16S profile, the path to the NCBI taxonomy and the output path,
downloads mapped genomes, creates an abundance table and all the further inputs which are needed downstream by the main pipeline.
If download is set to false, no genomes are downloaded and instead the reference genomes are expected to be in the out directory
The file name should then be out_path/taxID.fa.gz so it can be found
"""
#list of full genomes, input profile (CAMI format), taxonomy path, out directory
#returns number of genomes
def generate_input(args):
	genome_list, profile, tax_path, download, seed, no_samples, out_path, config, tax = read_args(args)

	tid,sci_name,ftp = read_genome_list(genome_list)
	
	profiles = transform_profile(profile,1,args.samples,tax) # might be multiple ones if biom file
	# probably 0.01 or something as threshold?
	
	downloaded, abundances, mapping, nr_samples = create_full_profiles(profiles, tid, ftp, tax, seed, download, out_path)

	return create_configs(nr_samples, out_path, config, abundances, downloaded, mapping)

