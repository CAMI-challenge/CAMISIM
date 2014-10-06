#!/usr/bin/env python

__author__ = 'hofmann'

import sys
import os
import argparse
import subprocess
from source.config import Config
from source.MetaTable import MetaTable


def ani_prediction_to_meta_table(metadata_table_in, metadata_table_out, input_querry_file, input_ref_file, mothur_otu_file, pool_size, ranks="family,genus,species"):
	executable = os.path.dirname(os.path.realpath(__file__)) + "/ani_prediction_to_meta_table.py"
	#print executable, "-i", metadata_table_in, "-o", metadata_table_out, "-cl", mothur_otu_file
	subprocess.call([executable, "-i", metadata_table_in, "-o", metadata_table_out, "-ir", input_ref_file, "-iq", input_querry_file, "-cl", mothur_otu_file, "-r", ranks, "-p", str(pool_size)])


def cluster_prediction_to_meta_table(metadata_table_in, metadata_table_out, mothur_cluster_file, threshold):
	executable = os.path.dirname(os.path.realpath(__file__)) + "/cluster_prediction_to_meta_table.py"
	#print executable, "-i", metadata_table_in, "-o", metadata_table_out, "-cl", mothur_otu_file
	subprocess.call([executable, "-i", metadata_table_in, "-o", metadata_table_out, "-cl", mothur_cluster_file, "-th", str(threshold)])


def detect_marker_genes(input_file, input_reference_file, output_folder, processors, is_ref_fasta_available=False):
	#folder_home = os.path.dirname(os.path.realpath(__file__))
	#folder_marker_genes = folder_home + "/marker_genes"
	executable = os.path.dirname(os.path.realpath(__file__)) + "/gather_16s.sh"
	if is_ref_fasta_available:
		#print executable, str(input_file), str(input_reference_file), output_folder, str(processors), "true"
		subprocess.call([executable, str(input_file), str(input_reference_file), output_folder, str(processors), "true"])
	else:
		#print executable, str(input_file), str(input_reference_file), output_folder, str(processors)
		subprocess.call([executable, str(input_file), str(input_reference_file), output_folder, str(processors)])


def cluster_marker_genes(output_folder, alignments_reference_db, threshold, distance_method, processors):
	suffix23s = "23S_rRNA.fna"
	suffix16s = "16S_rRNA.fna"
	suffix05s = "5S_rRNA.fna"

	executable = os.path.dirname(os.path.realpath(__file__)) + "/cluster.sh"

	input_file = output_folder + "/" + suffix16s
	if not os.path.isfile(input_file):
		print "Error -wd: File not found in working directory:", input_file
		return False
	#TODO: globaly declare such file name
	output_file = output_folder + "/" + suffix16s + ".otu.txt"
	#print "\n", executable, str(input_file), str(output_file), str(alignments_reference_db), str(threshold), str(distance_method), str(processors), "\n"
	subprocess.call([executable, str(input_file), str(output_file), str(alignments_reference_db), str(threshold), str(distance_method), str(processors)])

	#input_file = output_folder + "/" + suffix05s
	#if os.path.isfile(input_file):
	#	output_file = output_folder + "/" + suffix05s + ".otu.txt"
	#	#alignments_reference_db = "reference_db/nobackup/silva111/db/5S_bact+arch_dna.fna"
	#	print "./2_cluster.sh", str(input_file), str(output_file), str(alignments_reference_db), str(threshold), str(distance_method), str(processors)
	#	subprocess.call(['./2_cluster.sh', str(input_file), str(output_file), str(alignments_reference_db), str(threshold), str(distance_method), str(processors)])
	return True

def marker_gene_extraction(parser):
	"""The first step is to find and extract 16S marker gene sequences. The sequences are found using "hmmsearch" and extracted based on the given positions.
Two hmmer can currently be used. HMMER3.0 with a profile from 2010 and "rnammer" using HMMER2.0 with a profile from 2006.
A third one using HMMER3.0 is still to be evaluated.
So far it seems like rnammer provides better(more) results, but is also very slow.
input:
- file containing a list of fasta file paths, for example the genomes that someone wants to cluster.
- file containing a list of reference fasta file paths or alternatively, a fasta formated file containing the already extracted marker genes of the reference genomes.
- working directory where the results will be saved (intermediate files will be worked on in the designated /tmp folder)
- the number of processors that are available for parallel processing. The program "parallel" is used to process several genomes at the same time.
output:
- fasta formatted file containing the extracted marker genes of all genomes
"""
	args = parser.parse_args()

	input_file = args.input_file
	if input_file is None:
		print "Error -i: Please pass a file containing a list of genome paths: <id>\t<path>"
		return False

	if not os.path.isabs(input_file):
		input_file = os.getcwd() + "/" + input_file

	input_reference_fna_file = args.input_reference_fna_file
	input_reference_file = args.input_reference_file
	if input_reference_fna_file is None and input_reference_file is None:
		print "Error -ir(f): A reference file must be provided"
		return False
	#elif input_reference_fna_file is not None and input_reference_file is not None:
	#	print "Error -ir(f): too many arguments, -ir or -irf allowed, not both"
	#	return False

	if input_reference_file is not None:
		if not os.path.isabs(input_reference_file):
			input_reference_file = os.getcwd() + "/" + input_reference_file

	if input_reference_fna_file is not None:
		if not os.path.isabs(input_reference_fna_file):
			input_reference_fna_file = os.getcwd() + "/" + input_reference_fna_file

	working_directory = args.working_directory
	if working_directory is None:
		print "Error -wd: Please pass a working directory for results of the marker gene extraction"
		return False
	else:
		if not os.path.isabs(working_directory):
			working_directory = os.getcwd() + "/" + working_directory

	processors = args.processors
	if processors < 1:
		print "Error -p: only positive number of processors legal"
		return False

	if input_reference_fna_file is None:
		detect_marker_genes(input_file, input_reference_file, working_directory, processors)
	else:
		detect_marker_genes(input_file, input_reference_fna_file, working_directory, processors, True)
	return True


def gene_alignment_and_clustering(parser):
	"""The second step is to align 16S sequences and clustering them.
All is done using "mothur".
The alignment requires a high quality reference alignment (e.g. from SILVA) and
is done using the "gotoh" algorithm. (needleman, and blast are possible alternatives, but not implemented here)
Also using "mothur", empty or uninformative columns are removed.
When calculating distances (similar to DNADist) multi nucleotide gaps will be counted as one, gaps at the end are ignored.
To add even more references, the distances of the reference alignment will be merged with those of the working data.
These were precalculated and only the missing distances between the working data to the reference alignment need to be calculated.
The clustering will be done based on the final distances using the "Furthest neighbor" algorithm.
Mothur outputs cluster in steps up to a cutoff. The size of the steps can be chosen, by default 0.01 steps.
input:
- fasta formatted reference alignment (e.g. from SILVA)
- a threshold (between 0 and 1), that will be used for the clustering. 0.03 is default.
- working directory where the results will be saved and which contains the fasta formatted file with the extracted marker genes of all genomes
- the number of processors that are available for parallel processing. The program "mothur" can use several cores for the alignments and distance calculations.
output:
- a mothur formatted file containing the clusters, from unique up to the given threshold
"""
	args = parser.parse_args()

	working_directory = args.working_directory
	if working_directory is None:
		print "Error -wd: Please pass a working directory containing 16S genes for results of the clustering"
		return False
	else:
		if not os.path.isabs(working_directory):
			working_directory = os.getcwd() + "/" + working_directory

	processors = args.processors
	if processors < 1:
		print "Error -p: only positive number of processors legal"
		return False

	#alignments_reference_db = args.alignments_reference_db
	#if alignments_reference_db is None:
	#	print "Error -db: Please pass a file of markergenes in multi fasta format"
	#	return False

	distance_method = 2
	#if distance_method < 0 or distance_method > 2:
	#	print "Error -dist: Bad distance method, 0 phyml BioNJ, 1 raxml, 2 mothur"
	#	return False

	threshold = args.threshold

	folder_tools = os.path.dirname(os.path.realpath(__file__)) + "/tools"
	file_config = folder_tools + "/config.cfg"
	with open(file_config, 'r') as config_handler:
		config = Config(config_handler, 'PhyloPythiaS_Plus')
	#config = Config(folder_tools + "/config.cfg", 'PhyloPythiaS_Plus')
	folder_database = config.get('databaseFileSilver')
	if not os.path.isabs(folder_database):
		folder_database = os.path.normpath(folder_tools + "/" + folder_database)

	if folder_database is None:
		print "The taxonomy (databaseFile) is not specified."
		return False

	alignments_reference_db = None
	if os.path.isdir(folder_database):
		alignments_reference_db = os.path.join(os.path.normpath(folder_database), 'mothur_alignment_ref.fasta')
		#alignments_reference_db = os.path.join(os.path.normpath(folder_database), 'silva.archaea_bacteria.fasta')
		if not os.path.isfile(alignments_reference_db):
			print("The directory '{}' doesn't contain the reference file 'mothur_alignment_ref.fasta'".format(folder_database))
			return False

	if not os.path.isfile(alignments_reference_db):
		print("SILVA database file not found at:", alignments_reference_db)
		return False

	return cluster_marker_genes(working_directory, alignments_reference_db, threshold, distance_method, processors)


def classification_of_genomes_and_novelty_prediction(parser):
	"""As the third step, the unpublished genomes are classified based on the clusters they are found in.
Since clusters were made in 0.01 distance steps, the classification can be done using the smallest clusters first, using bigger ones if a classification can not be made.
If a marker gene of an unpublished genome is found in a cluster together with references, a common taxon that 90% of sequences agree with will be the predicted taxon.
The 90% is arbitrary chosen and is required because of taxonomic inconsistencies.
When a specific rank is checked for agreement, sequences with unknown classification on that rank are ignored.
TODO: check for taxonomic consitency on higher ranks for those!
Novelty prediction is based on the predicted taxon's rank. a high rank (phylum, order, class) with low distance can be a strong indicator for taxonomic inconsistencies.
But it could also be caused by sequences that are not fully classified, yet.
input:
- meta data table with a list of the genomes that are to be classified
- working directory where the results will be saved and which contains the mothur formatted file with the clusters
output:
- meta data table with a list of the genomes, with columns added that contain cluster based tax prediction, rank and novelty prediction
"""
	args = parser.parse_args()

	working_directory = args.working_directory
	if working_directory is None:
		print "Error -wd: Please pass a working directory containing the 16S genes clustering"
		return False
	else:
		if not os.path.isabs(working_directory):
			working_directory = os.getcwd() + "/" + working_directory

	metadata_table_out = args.metadata_table_out
	if metadata_table_out is None:
		print "Error -om: Please pass a output file path for the modified metatable"
		return False
	else:
		if not os.path.isabs(metadata_table_out):
			metadata_table_out = os.getcwd() + "/" + metadata_table_out

	metadata_table_in = args.metadata_table_in
	if metadata_table_in is None:
		metadata_table_in = metadata_table_out
		input_file = args.input_file
		if input_file is None or not os.path.isfile(input_file):
			print "Error -im: Please pass a metatable as input"
			return False
		else:
			if not os.path.isabs(input_file):
				input_file = os.getcwd() + "/" + input_file
		if not create_meta_table_from_fasta_path_file(input_file, metadata_table_in):
			print "Error -i: could not read file input"
			return False
	else:
		if not os.path.isabs(metadata_table_in):
			metadata_table_in = os.getcwd() + "/" + metadata_table_in

	#ncbitax_sqlite_db = args.ncbitax_sqlite_db
	#if ncbitax_sqlite_db is None:
	#	print "Error -tax: Please pass a path to a ncbi sqlite database"
	#	return False
	#else:
	#	if not os.path.isabs(ncbitax_sqlite_db):
	#		ncbitax_sqlite_db = os.getcwd() + "/" + ncbitax_sqlite_db

	#suffix16s = "16S_rRNA.fna"
	#TODO: declare this file name somewhere else, once! change name otu -> cluster
	mothur_cluster_file = working_directory + "/mothur_otu.txt"
	#otu_prediction_to_meta_table(metadata_table_in, metadata_table_out, mothur_otu_file, ncbitax_sqlite_db)

	cluster_prediction_to_meta_table(metadata_table_in, metadata_table_out, mothur_cluster_file, args.threshold)
	return True


def ani_of_genomes_and_novelty_prediction(parser):
	"""The fourth step is to calculate the average nucleotide identity.
To lessen the calculation burden, only genomes of sequences within the same cluster as an unpublished genome (marker gene) are compared.
In case only SILVA sequences are in a cluster, no comparison can be done.
The tool Mummer is used for the genome comparison, specifically nucmer.
Novelty predictions are made only for genomes with ani's better than 96%
ani > 96% -> same species
ani > 98% -> same strain
input:
- meta data table with a list of the genomes and data of the previous step, the output table path will be used as input.
- working directory which contains the mothur formatted file with the clusters (mothur_otu.txt)
- file containing a list of fasta file paths, the genomes that someone wants to cluster.
- file containing a list of reference fasta file paths.
output:
- meta data table with a list of the genomes, with columns added that contain ani based tax prediction, rank and novelty prediction
"""
	args = parser.parse_args()

	working_directory = args.working_directory
	if working_directory is None:
		print "Error -wd: Please pass a working directory containing the mothur formatted clusters file"
		return False

	if not os.path.isabs(working_directory):
		working_directory = os.getcwd() + "/" + working_directory

	metadata_table_out = args.metadata_table_out
	if metadata_table_out is None:
		print "Error -om: Please pass a output file path for the modified metatable"
		return False

	if not os.path.isabs(metadata_table_out):
		metadata_table_out = os.getcwd() + "/" + metadata_table_out
	metadata_table_in = metadata_table_out
	#metadata_table_in = args.metadata_table_in
	#if metadata_table_in is None:
	#	metadata_table_in = metadata_table_out
	#	if not os.path.isfile(metadata_table_out):
	#		print "Error -im: Please pass a metatable as input"
	#		return False
	#else:
	#	if not os.path.isabs(metadata_table_in):
	#		metadata_table_in = os.getcwd() + "/" + metadata_table_in

	input_querry_file = args.input_file
	if input_querry_file is None:
		print "Error -i: Please pass a file containing a list of genome paths: <id>\t<path>"
		return False

	if not os.path.isabs(input_querry_file):
		input_querry_file = os.getcwd() + "/" + input_querry_file

	input_reference_file = args.input_reference_file
	if input_reference_file is None:
		print "Error -ir: A file containing the file paths to the reference genomes is required"
		return False

	if input_reference_file is not None:
		if not os.path.isabs(input_reference_file):
			input_reference_file = os.getcwd() + "/" + input_reference_file

	pool_size = args.processors
	if pool_size < 1:
		print "Error -p: only positive number of processors legal"
		return False

	#TODO: declare this file name somewhere else, once! change name otu -> cluster
	mothur_cluster_file = working_directory + "/mothur_otu.txt"

	ani_prediction_to_meta_table(metadata_table_in, metadata_table_out, input_querry_file, input_reference_file, mothur_cluster_file, pool_size)
	return True


def create_meta_table_from_fasta_path_file(file_location_input, file_location_output):
	metadata_table = MetaTable(file_location_input, head=False)
	if metadata_table.number_of_rows == 0:
		return False
	id_column = metadata_table.get_column(0)
	#print len(id_column)
	#print len(metadata_table)
	metadata_table.clear()
	metadata_table.set_column("ID", id_column)
	metadata_table.write_meta_table_by_column(file_location_output)
	return True


def my_main():
	"""16S marker genes extraction, clustering, classification and novelty prediction

DESCRIPTION
Script for extracting 16S rRNA sequences from a set of genomes (new ones and reference genomes),
clustering of 16S sequences by first doing a multiple alignment, and then calculation the distance.
If the marker genes of the unpublished genomes clustered with one or more of the reference genomes, a taxonomic id is assign to them.
Finally the ani is calculated for those unpublished genomes, that got clustered with reference genomes.

The three stages in short:
1. Extraction of 16S sequences
2. Creation of a multiple alignment of 16S sequences, distance matrix calculation and clustering
3. Classification of genomes based on clustering and novelty prediction.
4. ANI calculation and novelty prediction.

input:
- file containing a list of fasta file paths, for example the genomes that someone wants to cluster.
- file containing a list of reference fasta file paths. (Step3) alternatively , a fasta formated file containing the already extracted marker genes of the reference genomes.
- a threshold (between 0 and 1), that will be used for the clustering. 0.03 is default.
- meta data table with a list of the genomes that are to be classified
- working directory where results of each stage will be saved
- the number of processors that are available for parallel processing.
output:
- meta data table with a list of the genomes, with columns added that contain tax predictions, average nucleotide identity and novelty predictions

TODO
- Evaluate third HMMER option
- confusion matrix of different cutoffs
- have all logfiles in the working directory
"""


	# example:
	usage = '''Example usage:
{0} \
-i id_to_path_unknown_genomes.txt \
-irf 16S_rRNA_refereces.fna \
-wd "my_output/" \
-th 0.10 \
-im metadata_table_test.csv \
-om metadata_table_test.out.csv \
-p 45 \
-s 0
'''.format(sys.argv[0])
	description = "taxonomic classify genomes by obtaining cluster based on 16S marker genes from genomes"
	parser = argparse.ArgumentParser(description=description)
	#parser = argparse.ArgumentParser(description=description, epilog=epilog)

	#parser.add_argument("-i", "--input_file", default=None, type=str,
	#					help="path to file containing marker genes in multi fasta format")

	#group = parser.add_mutually_exclusive_group()
	parser.add_argument("-irf", "--input_reference_fna_file", default=None, type=str,
						help="path to a fasta file containing the 16S marker genes of the reference genomes")
	parser.add_argument("-ir", "--input_reference_file", default=None, type=str,
						help="path to a file containing list of reference genomes file paths: <id>\\t<path>. NO COLUMN NAMES!")
	parser.add_argument("-i", "--input_file", default=None, type=str,
						help="path to a file containing tab separated list of genomes and their file paths: <id>\\t<path>. NO COLUMN NAMES!")

	parser.add_argument("-wd", "--working_directory", default=None, type=str,
						help="folder that will contain the log-files, found marker genes and also a file in mothur format containing the clustering")
	parser.add_argument("-th", "--threshold", default=0.03, type=float,
						help="threshold for clustering. Default: 0.03")
	#parser.add_argument("-db", "--alignments_reference_db", default=None, type=str,
	#					help="path to a 16S alignments reference database file (Silva)")

	parser.add_argument("-p", "--processors", default=2, type=int,
						help="number of processors to be used. Default: 2")
	parser.add_argument("-s", "--step", default=0, type=int,
						help='''available options: 0-4:
0 -> Full run through,
1 -> Marker gene extraction,
2 -> Gene alignment and clustering,
3 -> Classification of Genomes and novelty prediction
4 -> Average Nucleotide Identity calculation
Default: 0
''')

	parser.add_argument("-im", "--metadata_table_in", default=None, type=str,
						help="path to file containing tab separated list of genomes and their file path")
	parser.add_argument("-om", "--metadata_table_out", default=None, type=str,
						help="path to file containing tab separated list of genomes and their file path")
	#parser.add_argument("-tax", "--ncbitax_sqlite_db", default=None, type=str,
	#					help="path to file containing tab separated list of genomes and their file path")

	args = parser.parse_args()

	if args.step == 0 or args.step == 1:
		if not marker_gene_extraction(parser):
			#parser.print_help()
			sys.exit(1)

	if args.step == 0 or args.step == 2:
		if not gene_alignment_and_clustering(parser):
			#parser.print_help()
			sys.exit(1)

	if args.step == 0 or args.step == 3:
		if not classification_of_genomes_and_novelty_prediction(parser):
			#parser.print_help()
			sys.exit(1)

	if args.step == 0 or args.step == 4:
		if not ani_of_genomes_and_novelty_prediction(parser):
			#parser.print_help()
			sys.exit(1)

	sys.exit(0)

if __name__ == "__main__":
	my_main()
