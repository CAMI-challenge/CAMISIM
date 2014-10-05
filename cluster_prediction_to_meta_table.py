#!/usr/bin/env python

__author__ = 'hofmann'

from source.MothurOTU import MothurOTU
from source.MetaTable import MetaTable
from source.TaxonomyNcbi import TaxonomyNcbi
from source.config import Config
import sys
import os
import argparse


def my_main():
	"""Parsing of arguments"""
	# example:
	epilog = '''
	./tools/otu_prediction_to_meta_table.py
	-i metadata_table_v7.csv
	-o metadata_table_v7.out.csv
	-cl marker_genes/16S_rRNA.fna.otu.txt..otu.txt
	-db reference_db/nobackup/2014_02_01_ncbitax_sqlite.db
'''

	description = "Using otu clustering to make a prediction that is added to a metatable"
	parser = argparse.ArgumentParser(description=description)
	#parser = argparse.ArgumentParser(description=description, epilog=epilog)

	#parser.add_argument("-i", "--input_file", default=None, type=str,
	#					help="path to file containing marker genes in multi fasta format")
	parser.add_argument("-i", "--input_file", default=None, type=str,
						help="path to file containing tab separated list of genomes and their file path")
	parser.add_argument("-o", "--output_file", default=None, type=str,
						help="result folder that will contain OTUs in mothur format")
	parser.add_argument("-cl", "--cluster_file", default=None, type=str,
						help="path to file containing tab separated list of genomes and their file path")
	#parser.add_argument("-db", "--taxonomy_db", default=None, type=str,
	#					help="path to a taxonomy sqlite database file")
	parser.add_argument("-th", "--threshold", default=0.02, type=float,
						help="threshold for clustering. Default: 0.03")
	parser.add_argument("-minth", "--min_threshold", default=0.03, type=float,
						help="minimum threshold for clustering. Default: 0.02")
	args = parser.parse_args()

	if args.input_file is None:
		print "Error -i: Please pass a file containing a list of genome paths: <id>\t<path>"
		#parser.print_help()
		sys.exit(1)
	else:
		input_file = args.input_file
		if not os.path.isabs(input_file):
			input_file = os.getcwd() + "/" + input_file

	if args.output_file is None:
		print "Error -o: Please pass a file location for the results"
		#parser.print_help()
		sys.exit(1)
	else:
		output_file = args.output_file
		if not os.path.isabs(output_file):
			#output_file = os.path.realpath(__file__)+"/"+output_file
			output_file = os.getcwd() + "/" + output_file

	cluster_file = args.cluster_file
	if cluster_file is None or not os.path.isfile(cluster_file):
		print "Error -cl: Please pass a file containing a list of clusters in mothur format"
		#parser.print_help()
		sys.exit(1)
	if not os.path.isfile(cluster_file):
		print "Error -cl: File does not exist:", cluster_file
		#parser.print_help()
		sys.exit(1)

	if not os.path.isabs(cluster_file):
		cluster_file = os.getcwd() + "/" + cluster_file

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

	#print input_file
	#if args.step == 0 or args.step == 1:
	#	detect_marker_genes(input_file, output_folder, processors)
	#if args.step == 0 or args.step == 2:
	#	cluster_marker_genes(input_file, output_folder, alignments_reference_db, threshold, distance_method, processors)

	unpublished_genomes_id_column_name = "genome_ID"
	otu_id_column_name = "OTU_16S"
	otu_prediction_column_name = "OTU_NCBI_TAXONOMIC_PREDICTION"
	otu_prediction_sname_column_name = "OTU_SCIENTIFIC_NAME"
	otu_prediction_novelty_column_name = "OTU_NOVELTY_PREDICTION"
	otu_cutoff_column_name = "cutoff"

	print "reading metadata table file:", input_file
	metadata_table = MetaTable(input_file)
	metadata_table.remove_empty_columns()
	print "loading sqlite database:", taxonomy_db
	taxonomy = TaxonomyNcbi(taxonomy_db)
	print "reading mothur otu file:", cluster_file
	mothur_otus = MothurOTU(cluster_file, taxonomy)
	#print sorted(mothur_otus.otu_lists_by_cutoff_raw)
	#print mothur_otus.to_string("0.002")
	#sys.exit()

	otu_cutoff_column = metadata_table.get_column(otu_cutoff_column_name)
	if otu_cutoff_column is None:
		otu_cutoff_column = [''] * metadata_table.number_of_rows

	otu_id_column = metadata_table.get_column(otu_id_column_name)
	if otu_id_column is None:
		otu_id_column = [''] * metadata_table.number_of_rows

	otu_ncbi_column = metadata_table.get_column(otu_prediction_column_name)
	if otu_ncbi_column is None:
		otu_ncbi_column = [''] * metadata_table.number_of_rows

	otu_name_column = metadata_table.get_column(otu_prediction_sname_column_name)
	if otu_name_column is None:
		otu_name_column = [''] * metadata_table.number_of_rows

	otu_novelty_column = metadata_table.get_column(otu_prediction_novelty_column_name)
	if otu_novelty_column is None:
		otu_novelty_column = [''] * metadata_table.number_of_rows

	unpublished_genomes_id_column = metadata_table.get_column(unpublished_genomes_id_column_name)
	number_of_unpublished = len(unpublished_genomes_id_column)
	min_cutoff = str(args.min_threshold)
	all_done = False
	for cluster_cutoff in sorted(mothur_otus.lists_of_cluster_by_cutoff_raw):
		if all_done:
			break
		if cluster_cutoff == "unique":
			continue
		if float(cluster_cutoff) < float(min_cutoff):
			continue
		print "\n\n#cutoff", cluster_cutoff
		#print "\nget {0} ncbis of ids in otus using cutoff {1}".format(rank, otu_cutoff)
		mothur_otus.set_cutoff_raw(cluster_cutoff)
		all_done = True
		for row_index in range(0, number_of_unpublished):
			unpublished_genome_ids = unpublished_genomes_id_column[row_index]
			#metadata_table.get_entry_index(unknown_genomes_id_column_name, unknown_genomes_id)
			if unpublished_genome_ids != "" and otu_cutoff_column[row_index] == "" and mothur_otus.genome_id_exists(unpublished_genome_ids):
				all_done = False
				separator = ""
				#print "{0} of {1}: {2}".format(row_index + 1, number_of_unknown, unknown_genomes_id)
				list_of_cluster_id, list_of_cluster = mothur_otus.get_cluster_of(unpublished_genome_ids)
				if len(list_of_cluster_id) > 1:
					separator = ";"
				otu_ncbi = []
				otu_name = []
				otu_novelty = []
				otu_id = []
				for index in range(0, len(list_of_cluster_id)):
					cluster_id = list_of_cluster_id[index]
					cluster = list_of_cluster[index]
					if cluster is not None:
						ncbi_prediction, novelty = mothur_otus.get_cluster_ncbi_tax_prediction(cluster, unpublished_genomes_id_column, unpublished_genome_ids)
						if ncbi_prediction is not None and ncbi_prediction not in otu_ncbi:
							otu_cutoff_column[row_index] = str(cluster_cutoff)
							otu_ncbi.append(ncbi_prediction)
							otu_name.append(taxonomy.getScientificName(ncbi_prediction))
							otu_novelty.append("new_" + novelty)
							otu_id.append(str(ncbi_prediction) + "." + str(cluster_id))
				otu_ncbi_column[row_index] = separator.join(otu_ncbi)
				otu_name_column[row_index] = separator.join(otu_name)
				otu_novelty_column[row_index] =separator.join(otu_novelty)
				otu_id_column[row_index] = separator.join(otu_id)
						#else:
						#	print otu, "unknown species ncbi"
					#else:
						#otu_ncbi_column[row_index] = 1
						#print "no otu found for:", unknown_genomes_id

	metadata_table.set_column(otu_id_column_name, otu_id_column)
	metadata_table.set_column(otu_cutoff_column_name, otu_cutoff_column)
	metadata_table.set_column(otu_prediction_column_name, otu_ncbi_column)
	metadata_table.set_column(otu_prediction_sname_column_name, otu_name_column)
	metadata_table.set_column(otu_prediction_novelty_column_name, otu_novelty_column)
	metadata_table.write_meta_table_by_column(output_file)
	taxonomy.close()
	print "otu prediction finished"
	sys.exit(0)


if __name__ == "__main__":
	#help(start_process)
	my_main()
