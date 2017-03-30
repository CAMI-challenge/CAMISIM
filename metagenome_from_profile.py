#!/usr/bin/env python

from scripts.Validator.validator import Validator
from scripts.configfilehandler import ConfigFileHandler
import scripts.get_genomes as GG
import shutil
import os
import argparse

def parse_options():
	parser = argparse.ArgumentParser()
	
	#TODO default profile
	parser.add_argument("-p","--profile", default=None, type=str)

	# download the mapped full genomes?
	parser.add_argument("-dl", "--download-genomes", action='store_true', default=True)

	# out path
	parser.add_argument("-o", default=None, type=str)

	# temporary path
	parser.add_argument("-tmp", default=None, type=str)

	# file pointing to reference genomes TODO format
	parser.add_argument("-ref","--reference-genomes", default="tools/assembly_summary_complete_genomes.txt")

	# optional config file (out_path will get overwritten if it is set in config file)
	parser.add_argument("-c","--config",default="default_config.ini")

	parser.add_argument("-ncbi",default="tools/ncbi-taxonomy_20150130.zip")

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
	name = os.path.join(args.o,'') + "config.ini"
	with open(name,'wb') as nconf:
		nconf.write(new_config)
	return name

args = parse_options()
numg = GG.generate_input(args.reference_genomes,args.profile,args.ncbi,args.o,args.download_genomes)
c = create_config(args,numg)
os.system("./metagenomesimulation.py %s" % c)
