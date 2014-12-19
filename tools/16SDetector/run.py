#!/usr/bin/env python

"""
	Master script of the PPSplus.
"""

import sys
import os
import datetime
import argparse
import shutil
from algbioi.com import common
from algbioi.com import taxonomy_ncbi
from algbioi.com.config import Config
from algbioi.core.analysis_mg import MarkerGeneAnalysis
from algbioi.core.sequences import Sequences
from algbioi.core.analysis16s import RRNA16S
from algbioi.com import csv
from algbioi.ref import mask_db
from algbioi.core.pps import ppsOut2ppOut


def main():
	# external commands will be executed in Shell in Unix/Linux
	assert os.name == 'posix', str(
		'The pipeline runs only on "posix" systems (i.e. Unix/Linux compatible). ' + 'Your system is "' + os.name + '"')

	parser = argparse.ArgumentParser(
		description='''PhyloPythiaS Plus is an extension of PhyloPythiaS.''',
		epilog='''Read the user documentation for more details.''')

	parser.add_argument('-c', '--config', nargs=1, type=file, required=True, help='configuration file of the pipeline',
						metavar='config.cfg', dest='config')

	parser.add_argument('-i', '--input_fasta_file', nargs=1, default=None, required=False, help='path to fasta file',
						metavar='fasta_file.fn', dest='i')

	parser.add_argument('-out', '--output_folder', nargs=1, default=None, required=False, help='path to output folder',
						dest='out')

	parser.add_argument('-n', '--run-rrna16S', action='store_true',
						help='run hidden markov model and classify according to the 16S, 23S, and 5S marker genes',
						dest='n')

	parser.add_argument('-nn', '--runn-rrna16S', action='store_true',
						help='run hidden markov model searching for 16S, 23S, and 5S marker genes',
						dest='nn')

	parser.add_argument('-g', '--run-marker-gene-analysis', action='store_true',
						help='run hidden markov model and classify according to the "31" marker genes',
						dest='g')

	parser.add_argument("-hmmer", "--hmmer", default=3, type=int,
						help="'2': rnammer; '3': hmmsearch using hmmer 3.0")

	# push down predictions to more specific clades (sc)
	# build new OTUs (otu)
	# Taxator Blast LCA (jlcab), Taxator Last LCA (jlcal)
	# MlTreeMap (m); Amphora (ah); exclude some data from SSD (ex) place sequences;
	# (exSSD) exclude defined contigs from SSD; create training data'

	args = parser.parse_args()
	# read configuration
	config = Config(args.config[0], 'MarkerGeneExtraction')
	# pipeline directory
	if args.out is None:
		pipeline_dir = os.path.normpath(config.get('pipelineDir'))
	else:
		pipeline_dir = os.path.normpath(args.out[0])
	#print "output folder:", pipeline_dir
	if not os.path.isdir(pipeline_dir):
		print("Pipeline directory doesn't exist: ", pipeline_dir)
		return

	# create the following directories in the working directory if they don't exist
	working_dir = os.path.join(pipeline_dir, 'working')
	mg_working_dir = os.path.join(working_dir, 'mgWorking')
	output_dir = os.path.join(pipeline_dir, 'output')
	dir_array = [working_dir, output_dir, os.path.join(working_dir, 'projectDir'),
				 os.path.join(working_dir, 'sampleSpecificDir'), mg_working_dir]
	for dirPath in dir_array:
		if not os.path.isdir(dirPath):
			try:
				os.mkdir(dirPath)
			except OSError:
				print("Can't create directory", dirPath)
				return

	if args.i is None:
		input_fasta_file = os.path.normpath(config.get('inputFastaFile'))
	else:
		input_fasta_file = os.path.normpath(args.i[0])

	if (input_fasta_file is None) or (input_fasta_file is not None and (not os.path.isfile(input_fasta_file))):
		print("The input fasta file %s doesn't exist" % input_fasta_file)
		return
	if ('-' in input_fasta_file) or ('+' in input_fasta_file):
		print('The input fasta file path is not allowed to contain "+" or "-" characters ' +
			  'due to the "Mothur" software, given path:\n' + str(input_fasta_file))
		return

	input_fasta_scaffolds_file = None #config.get('inputFastaScaffoldsFile')
	if (input_fasta_scaffolds_file is not None) and (not os.path.isfile(input_fasta_scaffolds_file)):
		print("The given scaffolds fasta file doesn't exist: ", input_fasta_scaffolds_file)
		return

	scaffolds_to_contigs_map_file = None #config.get('scaffoldsToContigsMapFile')
	if (scaffolds_to_contigs_map_file is not None) and (not os.path.isfile(scaffolds_to_contigs_map_file)):
		print("The given scaffolds to contigs map file doesn't exist: ", scaffolds_to_contigs_map_file)
		return

	# create input id files (contigs, scaffolds, mappings)
	fasta_file_ids = common.createTagFilePath(working_dir, input_fasta_file, 'ids')
	seq_name_seq_id_file = common.createTagFilePath(working_dir, input_fasta_file, 'cToIds')

	fasta_file_scaffolds_ids = None
	scaff_name_scaff_id_file = None
	if input_fasta_scaffolds_file is not None:
		fasta_file_scaffolds_ids = common.createTagFilePath(working_dir, input_fasta_scaffolds_file, 'ids')
		scaff_name_scaff_id_file = common.createTagFilePath(working_dir, input_fasta_scaffolds_file, 'sToIds')

	scaffold_contig_map_ids_file = common.createTagFilePath(working_dir, input_fasta_file, 'mapSCIds')

	taxonomic_ranks = taxonomy_ncbi.TAXONOMIC_RANKS[1:]  # without root
	try:
		min_seq_len = int(config.get('minSeqLen'))
	except Exception as e:
		print("Can't parse configuration entry (minSeqLen), make sure that it's an integer number")
		raise e

	#####################################

	folder_tools = os.path.dirname(os.path.realpath(__file__)) + "/.."
	folder_database = config.get('databaseFile')
	if not os.path.isabs(folder_database):
		folder_database = os.path.normpath(folder_tools + "/" + folder_database)

	if folder_database is None:
		print("The taxonomy (databaseFile) is not specified.")
		return

	if not os.path.isdir(folder_database):
		print("Bad directory '{}'".format(folder_database))
		return

	database_file = os.path.join(os.path.normpath(folder_database), 'ncbitax_sqlite.db')
	if not os.path.isfile(database_file):
		print("The directory '%s' doesn't contain taxonomy file 'ncbitax_sqlite.db'")
		return

	if not os.path.isfile(database_file):
		print("The database file doesn't exist:", database_file)
		return

	# reference predictions
	reference_placement_file_out = None #config.get('referencePlacementFileOut')
	if reference_placement_file_out is not None and not os.path.isfile(reference_placement_file_out):
		print("The reference placement file doesn't exist: %s" % reference_placement_file_out)
		return

	# leaf clades contained in the reference
	if reference_placement_file_out is not None:
		try:
			refLeafCladesSet = set(map(int, csv.predToDict(reference_placement_file_out).values()))
		except Exception as e:
			print("Can't get clade taxon ids from the reference placement file: %s" % reference_placement_file_out)
			raise e
	else:
		refLeafCladesSet = None

	reference_placement_file_pp_out = None  # config.get('referencePlacementFilePPOut')
	if (reference_placement_file_pp_out is None) and (reference_placement_file_out is not None):
		# generate PP placement file
		try:
			reference_placement_file_pp_out = os.path.join(working_dir, os.path.basename(reference_placement_file_out) + '_.PP.out')
			ppsOut2ppOut(reference_placement_file_out, reference_placement_file_pp_out, taxonomic_ranks, database_file)
		except Exception as e:
			print("An error '%s' occurred while file '%s' has been generated."
				  % (e.message, reference_placement_file_pp_out))

	#####################################

	# generates working fasta files always when the configuration file is newer than particular files to be generated
	sequences = Sequences(input_fasta_file, input_fasta_scaffolds_file, scaffolds_to_contigs_map_file, taxonomic_ranks, min_seq_len)
	m_time_config = os.path.getmtime(config.getConfigFile())
	if ((not os.path.isfile(fasta_file_ids) or not os.path.isfile(seq_name_seq_id_file) or not os.path.isfile(scaffold_contig_map_ids_file)) or
			(m_time_config > min(os.path.getmtime(fasta_file_ids), os.path.getmtime(seq_name_seq_id_file), os.path.getmtime(scaffold_contig_map_ids_file)))):
		sequences.writeSequences(fasta_file_ids)
		print('Working contigs input fasta file created: %s' % fasta_file_ids)
		sequences.writeSeqNameSeqId(seq_name_seq_id_file)
		print('Ids mapping for the working contigs fasta file created: %s' % seq_name_seq_id_file)
		sequences.writeScaffoldContigMap(scaffold_contig_map_ids_file)
		print('Scaffolds -> contigs map ids file created: %s' % scaffold_contig_map_ids_file)
	# assert Common.seqFileCmp(inputFastaFile, fastaFileIds), 'The fasta IDs file contains different sequences!'
	if ((fasta_file_scaffolds_ids is not None) and ((not os.path.isfile(fasta_file_scaffolds_ids) or not os.path.isfile(scaff_name_scaff_id_file)) or
			(m_time_config > min(os.path.getmtime(fasta_file_scaffolds_ids), os.path.getmtime(scaff_name_scaff_id_file))))):
		sequences.writeScaffolds(fasta_file_scaffolds_ids)
		print('Working scaffolds input fasta file created: %s' % fasta_file_scaffolds_ids)  # (scaffold id -> sequence)
		sequences.writeScaffNameScaffId(scaff_name_scaff_id_file)
		print('Ids mapping for the working scaffolds fasta file created: %s' % scaff_name_scaff_id_file)

	# set up logging
	log_dir = os.path.join(output_dir, 'log')
	if not os.path.exists(log_dir):
		try:
			os.mkdir(log_dir)
		except OSError:
			print("Can't create logging directory:", log_dir)
			return

	# is it specified what to do?
	if not (args.n or args.nn or args.g or args.o or args.s):
		print('Choose what do you want to do!')
		print(parser.print_help())
		return

	# exclude mg from the reference ?
	s16_database = None  # os.path.normpath(config.get('s16Database'))
	exclude_ref_mg_rank = config.get('excludeRefMgRank')
	#s16_database = "/net/metagenomics/projects/cami_2014/02_otu_clustering/trunk/reference_db/nobackup/silva111" #config.get('s16Database')
	mg_database = None # os.path.normpath(config.get('mgDatabase'))

	if s16_database is not None and not os.path.isabs(s16_database):
		s16_database = os.path.normpath(folder_tools + "/" + s16_database)

	if mg_database is not None and not os.path.isabs(mg_database):
		mg_database = os.path.normpath(folder_tools + "/" + mg_database)

	#######################################################

	if (exclude_ref_mg_rank is not None) and (args.n or args.g):
		if not mask_db.isRankAllowed(exclude_ref_mg_rank):
			print("It's not allowed to exclude rank '%s' from the reference sequence database" % exclude_ref_mg_rank)
		else:
			mask_ref_mg_16s_dir = os.path.join(working_dir, str('ref_mg_16S_mask_' + exclude_ref_mg_rank))  # mask dir
			mask_ref_mg_dir = os.path.join(working_dir, str('ref_mg_mask_' + exclude_ref_mg_rank))  # mask dir
			mask_ref_ok = True
			if ((not os.path.isdir(mask_ref_mg_16s_dir)) or (not os.path.isdir(mask_ref_mg_dir)) or
					(min(os.path.getmtime(mask_ref_mg_16s_dir), os.path.getmtime(mask_ref_mg_dir)) < m_time_config)):
				if refLeafCladesSet is None:
					print("Can't exclude reference sequences since the reference file doesn't exist.")
					mask_ref_ok = False
				else:
					# remove the mask directories if they already exist
					if os.path.isdir(mask_ref_mg_16s_dir):
						shutil.rmtree(mask_ref_mg_16s_dir)
					if os.path.isdir(mask_ref_mg_dir):
						shutil.rmtree(mask_ref_mg_dir)
					# create the directories
					try:
						os.mkdir(mask_ref_mg_16s_dir)
						os.mkdir(os.path.join(mask_ref_mg_16s_dir, 'db'))
						os.mkdir(mask_ref_mg_dir)
						os.mkdir(os.path.join(mask_ref_mg_dir, 'db'))
					except OSError:
						print("Can't create one of the directories (%s, %s, 'db') for the masked reference marker gene "
							  "sequences." % (mask_ref_mg_16s_dir, mask_ref_mg_dir))
						mask_ref_ok = False
					else:
						mask_db.maskDb('mg', os.path.join(s16_database, 'db'), os.path.join(mask_ref_mg_16s_dir, 'db'),
									   exclude_ref_mg_rank, refLeafCladesSet, database_file)
						mask_db.maskDb('mg', os.path.join(mg_database, 'db'), os.path.join(mask_ref_mg_dir, 'db'),
									   exclude_ref_mg_rank, refLeafCladesSet, database_file)
						shutil.copy2(os.path.join(s16_database, 'content.csv'), mask_ref_mg_16s_dir)
						shutil.copy2(os.path.join(mg_database, 'content.csv'), mask_ref_mg_dir)
						shutil.copytree(os.path.join(mg_database, 'hmmAmphora'),
										os.path.join(mask_ref_mg_dir, 'hmmAmphora'))
			if mask_ref_ok:
				s16Database = mask_ref_mg_16s_dir
				mgDatabase = mask_ref_mg_dir
			else:
				print("Couldn't mask reference at rank '%s'." % exclude_ref_mg_rank)

	#######################################################

	# run 16S analysis
	if args.nn:
		hmmer = args.hmmer
		rrna = RRNA16S(config, s16_database, working_dir, folder_tools)
		print('\nrun Hidden Markov Model for 16S')
		rrna.runHMM(fasta_file_ids, outLog=get_log_file_name(log_dir, 'hmm16S'), hmmer=hmmer, moltypes="ssu", kingdoms="arc,bac")
		sys.exit(0)

	# run 16S analysis
	if args.n:
		hmmer = args.hmmer
		rrna = RRNA16S(config, s16_database, working_dir, folder_tools)
		print('\nrun Hidden Markov Model for (16S, 23S, 5S)')
		rrna.runHMM(fasta_file_ids, outLog=get_log_file_name(log_dir, 'hmm16S'), hmmer=hmmer, moltypes="lsu,ssu,tsu", kingdoms="arc,bac")
		print('\nrun (16S, 23S, 5S) classification')
		rrna.classify16S(fasta_file_ids, outLog=get_log_file_name(log_dir, 'mothurClassify16S'))
		rrna.classify23S(fasta_file_ids, outLog=get_log_file_name(log_dir, 'mothurClassify23S'))
		rrna.classify5S(fasta_file_ids, outLog=get_log_file_name(log_dir, 'mothurClassify5S'))

	# marker gene analysis (Amphora)
	if args.g:
		mg = MarkerGeneAnalysis(config, mg_database, working_dir, mg_working_dir)
		print 'run 31 marker gene analysis'
		mg.runMarkerGeneAnalysis(fasta_file_ids, outLog=get_log_file_name(log_dir, 'amphoraMG'))

	###############################################

	# remove non-DNA characters in the sequences of contigs and scaffolds
	# this is important before the sample specific data are created and before PPS predicts the data
	sequences.setRemoveNonDna(True)
	# the working fasta files of scaffolds and contigs are regenerated, so that they
	# wouldn't contain non-DNA characters
	sequences.writeSequences(fasta_file_ids)
	print('Working contigs input fasta file without non-DNA chars created: %s' % fasta_file_ids)
	if fasta_file_scaffolds_ids is not None:
		sequences.writeScaffolds(fasta_file_scaffolds_ids)
		print('Working scaffolds input fasta file without non-DNA chars created: %s' % fasta_file_scaffolds_ids)


def get_log_file_name(log_dir, description):
	now = datetime.datetime.now()
	time_stamp = str(str(now.year) + str(now.month) + str(now.day) + '_' + str(now.hour) + str(now.minute) + str(now.second))
	file_name = str(description) + '_' + time_stamp + '.txt'
	return os.path.normpath(os.path.join(log_dir, file_name))


if __name__ == "__main__":
	main()
