#!/usr/bin/env python

__author__ = 'hofmann'

from source.MothurOTU import MothurOTU
from source.MetaTable import MetaTable
from source.TaxonomyNcbi import TaxonomyNcbi
from source.config import Config
from source.ANIm import ANIm
import sys
import os
import argparse
import logging
import logging.handlers

#original prototype:
#   http://armchairbiology.blogspot.de/2013/11/ani-are-you-okay-are-you-okay-ani.html
#   (c) The James Hutton Institute 2013
#   Author: Leighton Pritchard
#
#   Contact:
#   leighton.pritchard@hutton.ac.uk
#   GNU General Public License


def get_logger(verbose=False, logfile=None):
	# We set up logging, and modify loglevel according to whether we need
	# verbosity or not
	# err_handler points to sys.stderr
	# err_handler_file points to a logfile, if named
	logger = logging.getLogger('calculate_ani.py')
	logger.setLevel(logging.DEBUG)
	err_handler = logging.StreamHandler(sys.stderr)
	err_formatter = logging.Formatter('%(levelname)s: %(message)s')
	err_handler.setFormatter(err_formatter)
	if logfile is not None:
		try:
			log_stream = open(logfile, 'w')
			err_handler_file = logging.StreamHandler(log_stream)
			err_handler_file.setFormatter(err_formatter)
			err_handler_file.setLevel(logging.INFO)
			logger.addHandler(err_handler_file)
		except:
			logger.error("Could not open %s for logging" % logfile)
			sys.exit(1)
	if verbose:
		err_handler.setLevel(logging.INFO)
	else:
		err_handler.setLevel(logging.WARNING)
	logger.addHandler(err_handler)
	return logger


def my_main():
	"""Parsing of arguments"""
	# example:
	epilog = '''
'''

	description = "Using average nucleotide identity to make novelty predictionw that are added to a metatable"
	parser = argparse.ArgumentParser(description=description)
	#parser = argparse.ArgumentParser(description=description, epilog=epilog)

	parser.add_argument("-i", "--metadata_table_in", default=None, type=str,
						help="")
	parser.add_argument("-o", "--metadata_table_out", default=None, type=str,
						help="")
	parser.add_argument("-iq", "--input_querry_file", default=None, type=str,
						help="path to file containing tab separated list of genomes and their file path")
	parser.add_argument("-ir", "--input_ref_file", default=None, type=str,
						help="path to file containing tab separated list of genomes and their file path")
	parser.add_argument("-cl", "--mothur_otu_file", default=None, type=str,
						help="")
	parser.add_argument("-p", "--pool_size", default=10, type=int,
						help="number of processors to be used. Default: 2")
	parser.add_argument("-r", "--ranks", default="family,genus,species", type=str,
						help="legal ranks, default: family,genus,species")
	args = parser.parse_args()

	metadata_table_in = args.metadata_table_in
	if metadata_table_in is None:
		print "Error -i: "
		#parser.print_help()
		sys.exit(1)

	metadata_table_out = args.metadata_table_out
	if metadata_table_out is None:
		print "Error -o: "
		#parser.print_help()
		sys.exit(1)

	input_querry_file = args.input_querry_file
	if input_querry_file is None:
		print "Error -iq: "
		sys.exit(1)

	input_ref_file = args.input_ref_file
	if input_ref_file is None:
		print "Error -ir: "
		sys.exit(1)

	mothur_otu_file = args.mothur_otu_file
	if mothur_otu_file is None:
		print "Error -cl: Please pass a file containing a list of clusters in mothur format"
		#parser.print_help()
		sys.exit(1)

	pool_size = args.pool_size
	if pool_size < 1:
		print "Error -p: only positive number of processors legal"
		sys.exit(1)

	if not os.path.isabs(metadata_table_in):
		metadata_table_in = os.getcwd() + "/" + metadata_table_in

	if not os.path.isabs(metadata_table_out):
		#output_file = os.path.realpath(__file__)+"/"+output_file
		metadata_table_out = os.getcwd() + "/" + metadata_table_out

	if not os.path.isabs(input_querry_file):
		input_querry_file = os.getcwd() + "/" + input_querry_file

	if not os.path.isabs(input_ref_file):
		input_ref_file = os.getcwd() + "/" + input_ref_file

	if not os.path.isabs(mothur_otu_file):
		mothur_otu_file = os.getcwd() + "/" + mothur_otu_file

	folder_tools = os.path.dirname(os.path.realpath(__file__)) + "/tools"
	file_config = folder_tools + "/config.cfg"
	#print file_config
	#config = None
	with open(file_config, 'r') as config_handler:
		config = Config(config_handler, 'PhyloPythiaS_Plus')
	folder_database = config.get('databaseFile')
	if not os.path.isabs(folder_database):
		folder_database = os.path.normpath(folder_tools + "/" + folder_database)

	if folder_database is None:
		print("The taxonomy (databaseFile) is not specified.")
		sys.exit(1)

	taxonomy_db = None
	if os.path.isdir(folder_database):
		taxonomy_db = os.path.join(os.path.normpath(folder_database), 'ncbitax_sqlite.db')
		if not os.path.isfile(taxonomy_db):
			print("The directory '%s' doesn't contain sqllite taxonomy database file 'ncbitax_sqlite.db'")
			sys.exit(1)

	if not os.path.isfile(taxonomy_db):
		print("The sqllite taxonomy database file doesn't exist:", taxonomy_db)
		sys.exit(1)

	ranks = args.ranks.strip().split(',')

	unpublished_genomes_id_column_name = "genome_ID"
	otu_cutoff_column_name = "cutoff"
	ani_prediction_novelty_column_name = "ANI_NOVELTY_PREDICTION"
	ani_prediction_column_name = "ANI_NCBI_TAXONOMIC_PREDICTION"
	ani_scientific_name_column_name = "ANI_SCIENTIFIC_NAME"
	ani_column_name = "ANI"

	print "reading metadata table file:", metadata_table_in
	metadata_table = MetaTable(metadata_table_in)
	print "loading sqlite database:", taxonomy_db
	taxonomy = TaxonomyNcbi(taxonomy_db)
	print "reading mothur otu file:", mothur_otu_file
	mothur_otus = MothurOTU(mothur_otu_file, taxonomy)
	#print sorted(mothur_otus.otu_lists_by_cutoff_raw)
	#print mothur_otus.to_string("0.002")
	#sys.exit()

	ani_scientific_name_column = metadata_table.get_column(ani_scientific_name_column_name)
	if ani_scientific_name_column is None:
		ani_scientific_name_column = [''] * metadata_table.number_of_rows

	otu_cutoff_column = metadata_table.get_column(otu_cutoff_column_name)
	if otu_cutoff_column is None:
		otu_cutoff_column = [''] * metadata_table.number_of_rows

	ani_prediction_novelty_column = metadata_table.get_column(ani_prediction_novelty_column_name)
	if ani_prediction_novelty_column is None:
		ani_prediction_novelty_column = [''] * metadata_table.number_of_rows

	ani_prediction_column = metadata_table.get_column(ani_prediction_column_name)
	if ani_prediction_column is None:
		ani_prediction_column = [''] * metadata_table.number_of_rows

	ani_column = metadata_table.get_column(ani_column_name)
	if ani_column is None:
		ani_column = [''] * metadata_table.number_of_rows

	unpublished_genome_ids_column = metadata_table.get_column(unpublished_genomes_id_column_name)

	logger = get_logger(True)
	nucmer_exe = "importpackage mummer;nucmer"
	with ANIm(input_querry_file, input_ref_file, None, nucmer_exe=nucmer_exe, logger=logger, pool_size=pool_size) as ani_calculator:
		otus = mothur_otus.get_clusters(otu_cutoff_column, unpublished_genome_ids_column)
		#logger.info("OTUS: {}".format(otus))
		#sys.exit(0)
		for unknown_genomes_id in unpublished_genome_ids_column:
			if otus[unknown_genomes_id] is None:
				continue
			ani_calculator.add_nucmer_cmd_lines(otus[unknown_genomes_id], [unknown_genomes_id])

		min_lengths, min_sim_errors, min_perc_ids, min_perc_aln, ncbi, ani_ish = ani_calculator.calculate_minimum_anim()

		number_of_unknown = len(unpublished_genome_ids_column)
		for row_index in range(0, number_of_unknown):
			unknown_genomes_id = unpublished_genome_ids_column[row_index]
			if unknown_genomes_id in ani_ish and ani_ish[unknown_genomes_id] > 0:
				if float(ani_ish[unknown_genomes_id]) > 0.98:
					ani_prediction_novelty_column[row_index] = "same_strain"
				elif float(ani_ish[unknown_genomes_id]) > 0.96:
					ani_prediction_novelty_column[row_index] = "same_species"
				ani_prediction_column[row_index] = ncbi[unknown_genomes_id]
				ani_column[row_index] = str(ani_ish[unknown_genomes_id])
				science_name = taxonomy.getScientificName(ncbi[unknown_genomes_id])
				if science_name is not None:
					ani_scientific_name_column[row_index] = science_name

	taxonomy.close()
	metadata_table.set_column(ani_column_name, ani_column)
	metadata_table.set_column(ani_prediction_novelty_column_name, ani_prediction_novelty_column)
	metadata_table.set_column(ani_prediction_column_name, ani_prediction_column)
	metadata_table.set_column(ani_scientific_name_column_name, ani_scientific_name_column)
	metadata_table.write_meta_table_by_column(metadata_table_out)
	print "ani prediction finished"
	sys.exit(0)


if __name__ == "__main__":
	my_main()
