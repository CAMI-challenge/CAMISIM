import sys
import os
from ftplib import FTP 
import gzip
import random
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy
from scripts.Validator.validator import Validator
from scripts.loggingwrapper import LoggingWrapper as logger

"""
Given a 16S-profile (currently only in CAMI format), downloads all closest relative genomes and creates abundances
"""

_ranks=['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
# strain inclusion?

"""
original code in the profiling-evaluation-biobox, reads a file in cami profiling format and extracts relevant information (taxids/tax path/relative abundance/genome rank in the taxonomy)
adpated from the profile evaluation biobox, extendeded by the following: We only check for the species tax ids, original genomes for the higher ranks will be checked later on
"""
# TODO generic read profile function, s.t. biom/QIIME format can also be used
def read_profile(file_path, epsilon):
	assert isinstance(file_path, basestring)
	if isinstance(file_path, str) and not os.path.isfile(file_path):
		logger.error("16S profile not found in: %s" % file_path)
		raise Exception("File not found")
	assert epsilon is None or isinstance(epsilon, (float, int, long))
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
	return tax_ids, tax_path, weights, ranks

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
		logger.error("Reference genome list not found in: %s" % file_path)
		raise Exception("File not found")
	tax_ids = list()
	sci_name = list()
	ftp_address = {}
	with open(file_path,'r') as full_genomes:
		for line in full_genomes:
			if len(line) == 0:
				continue
			temp = line.split('\t')
			if len(temp) < 3: # if the line split by tabs does not contain the three elements something is wrong
				continue
			tax_ids.append(temp[0])
			sci_name.append(temp[1])
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
	for rank in _ranks[1:]:
		ids_per_rank.append(dict()) # create empty dicts
	for tax_id in tax_ids:
		try:
			lineage = tax.get_lineage_of_legal_ranks(tax_id, ranks = _ranks)
		except ValueError:
			continue # this means, the taxid was not found in the reference (reference mismatch)
		i = 1
		for rank_tax_id in lineage[1:]: # for the path up to the source/superkingdom
			if rank_tax_id in ids_per_rank[i]:# add the species id to the map, if the higher rank tax id is present
				ids_per_rank[i][rank_tax_id].append(tax_id) 
			else:
				ids_per_rank[i].update({rank_tax_id : [tax_id]})
			i = i + 1 # continue with the next rank
	return ids_per_rank

""" 
Given the list of available full genomes and the species extracted from the profile, finds a mapping of profile genomes to available full genomes sequences
First, extends the list of reference tax ids, so the least common ancestors in the phylogeny can be found.
For all species tax ids in the profile, it is checked whether that tax id is a tax id of a known full genome.
If not, than the next rank is checked, i.e. genus. If there is a full genome with the same genus tax id like our genus tax id,
than one of these genomes is chosen as the "closest related" genome. 
"""
def map_to_full_genomes(ref_tax_ids, profile_tax_ids, tax, seed):
	random.seed(seed)
	gen_map = {}
	extended_genome_list = extend_genome_list(ref_tax_ids,tax)
	to_download = dict() # ncbi id of genomes to download
	for taxid in profile_tax_ids[0]:
		found_genome = False
		if taxid in extended_genome_list[0]: # a full genome with exact ncbi id is present
			to_download.update(taxid = taxid) # if profile contains strains this might cause overwrites TODO
			found_genome = True
		else: # the exact genome is not present, go up the ranks
			try:
				lineage = tax.get_lineage_of_legal_ranks(taxid,ranks = _ranks)
			except ValueError: #tax ID was not found in reference data base
				logger.warning("Genome %s not found in reference, maybe your reference is deprecated?" % taxid)
				continue
			i = 0
			for higher_taxid in lineage: # rank is a number corresponding to the ranks defined in _ranks with species being the lowers (0)
				if higher_taxid in extended_genome_list[i]:
					species_id = extended_genome_list[i][higher_taxid]
					to_download.update({taxid : species_id[random.randint(0,len(species_id) - 1)]}) #randomly select one of the mapped genomes TODO 
					found_genome = True
					break #TODO add rank for debugging purposes (_ranks[i])
				i = i + 1 # go to the next rank
		if not found_genome: # This should not happen, if tax id is in reference! TODO
			logger.warning("No genome corresponding to ID %s found, omitted. Maybe your reference is deprecated?" % taxid)
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
def download_genomes(list_of_genomes, ftp_list, out_path):
	metadata = dict() # create the metadata table (pathes to genomes)
	ftp = FTP('ftp.ncbi.nlm.nih.gov') 
	ftp.login() # anonymous login
	for elem in list_of_genomes:
		path = ftp_list[list_of_genomes[elem]]
		split_path = path.split('/')
		cwd = "/" + "/".join(split_path[3:]).rstrip() # get /address/to/genome
		gen_name = split_path[-1].rstrip() # genome name is last in address
		to_dl = gen_name + "_genomic.fna.gz"
		out_name = out_path + list_of_genomes[elem] + ".fa"  # out name is the ncbi id of the downloaded genome
		out_name_gz = out_name + ".gz"
		metadata.update({list_of_genomes[elem]:out_name})
		ftp.cwd(cwd)
		ftp.retrbinary("RETR %s" % to_dl, open(out_name_gz,'wb').write)
		gf = gzip.open(out_name_gz) #TODO we have to "uncompress" the files
		outF = open(out_name,'wb')
		outF.write(gf.read())
		gf.close()
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
		ab = int(ab*100) # just so we get nice numbers
		abundance.update({list_of_genomes[elem]:ab})
	return abundance

"""
Given the reference genomes' sequences (path), an 16S profile, the path to the NCBI taxonomy and the output path,
downloads mapped genomes, creates an abundance table and all the further inputs which are needed downstream by the main pipeline.
If download is set to false, no genomes are downloaded and instead the reference genomes are expected to be in the out directory
The file name should then be out_path/taxID.fa.gz so it can be found
"""
#list of full genomes, input profile (CAMI format), taxonomy path, out directory
#returns number of genomes
def generate_input(genome_list,profile,tax_path,out_path,download):
	out_path = os.path.join(out_path,'') #so we are sure it is a directory
	tid,n,ftp = read_genome_list(genome_list)
	profile = read_profile(profile,1) # probably 0.01 or something as threshold?
	tax = NcbiTaxonomy(tax_path)
	to_dl = map_to_full_genomes(tid,profile,tax,None) #TODO add seed
	if not download:
		downloaded = []
		for gen in to_dl:
			downloaded.append("%s\t%s" % (to_dl[gen],out_path + to_dl[gen] + "fa.gz"))
	else:
		downloaded = download_genomes(to_dl,ftp,out_path)
	ab = create_abundance_table(to_dl,profile)
	with open(out_path + "abundance.tsv",'wb') as abundance:
		for e in ab:
			abundance.write("%s\t%s\n" % (e,ab[e]))
	with open(out_path + "genome_to_id.tsv",'wb') as gpath:
		for e in downloaded:
			gpath.write("%s\t%s\n" % (e,downloaded[e]))
	with open(out_path + "metadata.tsv",'wb') as metadata:
		i = 0 #for OTU assignment (every species gets its own OTU here) TODO
		for gen in to_dl:
			metadata.write("%s\t%s\t%s\t%s\n" % (gen,i,to_dl[gen],"new_strain"))
			i = i + 1
	return len(downloaded)

