#!/usr/bin/env python

from scripts.Validator.validator import Validator
from scripts.configfilehandler import ConfigFileHandler
import scripts.get_genomes as GG
import shutil
import os
import argparse

def parse_options():
	parser = argparse.ArgumentParser()
	
	helptext="16S profile to create metagenome from. Can either be CAMI-format or biom-format."
	#TODO default profile
	parser.add_argument("-p","--profile", default=None, type=str,help=helptext)

	helptext="Whether the related genomes are supposed to be downloaded."
	# download the mapped full genomes?
	parser.add_argument("-dl", "--download-genomes", action='store_true', default=True, help=helptext)

	helptext="Output directory, make sure this directory exists!"
	# out path
	parser.add_argument("-o", default=None, type=str, help=helptext)
	
	helptext="Path where temporary files are stored (get deleted after pipeline is finished)"
	# temporary path
	parser.add_argument("-tmp", default=None, type=str, help=helptext)

	helptext="File pointing to reference genomes of the format: NCBI id\tScientific name\tNCBI ftp address of full genome"
	# file pointing to reference genomes TODO format
	parser.add_argument("-ref","--reference-genomes", default="tools/assembly_summary_complete_genomes.txt", help=helptext)

	helptext="Path to config file. Careful when setting \"metadata\", \"id_to_genome_file\", \"distribution_file_paths\"(they will be set by the pipeline) and the out path differently from the command line out path"
	# optional config file (out_path will get overwritten if it is set in config file)
	parser.add_argument("-c","--config",default="default_config.ini",help=helptext)

	helptext="Path to the NCBI taxdump for finding corresponding reference genomes"
	parser.add_argument("--ncbi",default="tools/ncbi-taxonomy_20150130.zip",help=helptext)
	
	helptext="Seed for the random generator"
	parser.add_argument("--seed",default=None,help=helptext)

	args = parser.parse_args()
	
	return args

def create_config(args,numg):
	with open(args.config) as config:
		new_config = config.read()
	
	new_config+="metadata=%s\n" % (os.path.join(args.o,'') + "metadata.tsv")
	new_config+="id_to_genome_file=%s\n" % (os.path.join(args.o,'') + "genome_to_id.tsv")
	new_config+="distribution_file_paths=%s\n" % (os.path.join(args.o,'') + "abundance.tsv")
	new_config+="output_directory=%s\n" % (os.path.join(args.o,''))
	new_config+="genomes_total=%s\n" % numg
	if args.seed is not None:
		new_config+="seed=%s\n" % args.seed
	name = os.path.join(args.o,'') + "config.ini"
	with open(name,'wb') as nconf:
		nconf.write(new_config)
	return name

args = parse_options()
numg = GG.generate_input(args.reference_genomes,args.profile,args.ncbi,args.o,args.download_genomes,args.seed)
c = create_config(args,numg)
os.system("./metagenomesimulation.py %s" % c)
