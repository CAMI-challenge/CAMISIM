#!/usr/bin/env python

from scripts.Validator.validator import Validator
from scripts.configfilehandler import ConfigFileHandler
try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser  # ver. < 3.0
import scripts.get_genomes as GG
import shutil
import os
import argparse

def parse_options():
	parser = argparse.ArgumentParser()
	
	helptext="16S profile to create metagenome from. Can either be CAMI-format or biom-format."
	#TODO default profile
	parser.add_argument("-p","--profile", default=None, type=str, help=helptext)

	helptext="Number of samples to be generated. If nothing is given, this defaults to 1 (CAMI format) or the number of samples present in the biom file. If a specific number is given, the samples are simulated using the first sample of the biom file"
	#TODO is this the wanted behaviour?
	parser.add_argument("-s","--samples", default=None, type=int, help=helptext)

	helptext="Whether the related genomes are supposed to be downloaded."
	# download the mapped full genomes?
	parser.add_argument("-no-dl", "--dont-download-genomes", action='store_false', default=True, help=helptext)

	helptext="Output directory, make sure this directory exists!"
	# out path
	parser.add_argument("-o", default=None, type=str, help=helptext,metavar="OUT PATH")
	
	helptext="Path where temporary files are stored (get deleted after pipeline is finished)"
	# temporary path
	parser.add_argument("-tmp", default=None, type=str, help=helptext)

	helptext="File pointing to reference genomes of the format: NCBI id\tScientific name\tNCBI ftp address of full genome"
	# file pointing to reference genomes TODO format
	parser.add_argument("-ref","--reference-genomes", default="tools/assembly_summary_complete_genomes.txt", help=helptext)

	helptext="Path to config file. Careful when setting \"metadata\", \"id_to_genome_file\", \"distribution_file_paths\"(they will be set by the pipeline) and the out path differently from the command line out path"
	# optional config file (out_path will get overwritten if it is set in config file)
	parser.add_argument("-c","--config",default="default_config.ini",help=helptext,metavar="CONFIG FILE")

	helptext="Path to the NCBI taxdump for finding corresponding reference genomes"
	parser.add_argument("--ncbi",default="tools/ncbi-taxonomy_20170222.tar.gz",help=helptext)
	
	helptext="Seed for the random generator"
	parser.add_argument("--seed",default=None,help=helptext)

	args = parser.parse_args()
	
	return args

def create_config(args,cfg,numg):
	config = ConfigParser()
	config.read(cfg)
	
	config.set('Main', 'output_directory', os.path.join(args.o,''))

	if args.seed is not None:
		config.set('Main', "seed", args.seed)
	name = os.path.join(args.o,'') + "config.ini"
	with open(name,'wb') as cfg_path:
		config.write(cfg_path)
	return name

args = parse_options()
numg,config = GG.generate_input(args) # total number of genomes and path to updated config
c = create_config(args,config,numg)
os.system("./metagenomesimulation.py %s" % c)
