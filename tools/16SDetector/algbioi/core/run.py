#!/usr/bin/env python

"""
	Master script of the PPSplus.
"""
import os
import datetime
import argparse
from algbioi.com import common
from algbioi.com import taxonomy_ncbi
from algbioi.com.config import Config
from algbioi.core.analysis_mg_ import MarkerGeneAnalysis
from algbioi.core.sequences import Sequences
from algbioi.core.analysis16s import RRNA16S

#import subprocess
#import traceback
#import shutil
#from algbioi.com import csv
#from algbioi.com import fasta as fas
#from algbioi.core import ssd_eval
#from algbioi.core.pps import generateCladesForGeneralModel
#from algbioi.core.training_data import PPSInput
#from algbioi.core.cluster import MGCluster
#from algbioi.core.taxonomy import Taxonomy
#from algbioi.eval import consistency
#from algbioi.eval import accuracy
#from algbioi.eval import confusion_matrix
#from algbioi.ref import mask_db
#from algbioi.core.pps import readPPSOutput
#from algbioi.core.pps import computeTrainingAccuracy
#from algbioi.misc import out_proc
#from algbioi.core.pps import ppsOut2ppOut

# paths on hera/gaia:
# export PATH=/net/programs/Debian-6.0.3-x86_64/python-2.7joh/bin:$PATH
# export PYTHONPATH=/net/metagenomics/projects/PPSmg/scripts/scriptsR30
#
from numpy.distutils.system_info import system_info


def main():
	# external commands will be executed in Shell in Unix/Linux
	assert os.name == 'posix', str(
		'The pipeline runs only on "posix" systems (i.e. Unix/Linux compatible). ' + 'Your system is "' + os.name + '"')

	parser = argparse.ArgumentParser(
		description='''PhyloPythiaS Plus is an extension of PhyloPythiaS.''',
		epilog='''Read the user documentation for more details.''')

	parser.add_argument('-c', '--config', nargs=1, type=file, required=True, help='configuration file of the pipeline',
						metavar='config.cfg', dest='config')

	parser.add_argument('-n', '--run-rrna16S', action='store_true',
						help='run hidden markov model and classify according to the 16S, 23S, and 5S marker genes',
						dest='n')

	parser.add_argument('-g', '--run-marker-gene-analysis', action='store_true',
						help='run hidden markov model and classify according to the "31" marker genes',
						dest='g')

	#parser.add_argument('-o', '--process-preprocessed-output', action='store', nargs='+',
	#					help='process output from the 16S rRNA analysis (s16), marker gene "31" analysis (mg); '
	#						 'to build the general model (general)',
	#					choices=["s16", "mg", "general"],
	#					dest='o')
	# push down predictions to more specific clades (sc)
	# build new OTUs (otu)
	# Taxator Blast LCA (jlcab), Taxator Last LCA (jlcal)
	# MlTreeMap (m); Amphora (ah); exclude some data from SSD (ex) place sequences;
	# (exSSD) exclude defined contigs from SSD; create training data'
	args = parser.parse_args()
	# read configuration
	config = Config(args.config[0], 'PhyloPythiaS_Plus')

	# pipeline directory
	pipeline_dir = os.path.normpath(config.get('pipelineDir'))
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

	# read input fasta: contigs, scaffolds, mapping
	input_fasta_file = os.path.normpath(config.get('inputFastaFile'))
	if (input_fasta_file is None) or (input_fasta_file is not None and (not os.path.isfile(input_fasta_file))):
		print("The input fasta file %s doesn't exist" % input_fasta_file)
		return
	if ('-' in input_fasta_file) or ('+' in input_fasta_file):
		print('The input fasta file path is not allowed to contain "+" or "-" characters ' +
			  'due to the "Mothur" software, given path:\n' + str(input_fasta_file))
		return

	input_fasta_scaffolds_file = config.get('inputFastaScaffoldsFile')
	if (input_fasta_scaffolds_file is not None) and (not os.path.isfile(input_fasta_scaffolds_file)):
		print("The given scaffolds fasta file doesn't exist: ", input_fasta_scaffolds_file)
		return

	scaffolds_to_contigs_map_file = config.get('scaffoldsToContigsMapFile')
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
	if not (args.n or args.g or args.o or args.t or args.p or args.r or args.s or args.a):
		print('Choose what do you want to do!')
		print(parser.print_help())
		return

	# exclude mg from the reference ?
	s16_database = os.path.normpath(config.get('s16Database'))
	mg_database = os.path.normpath(config.get('mgDatabase'))

	# run 16S analysis
	if args.n:
		rrna = RRNA16S(config, s16_database, working_dir)
		print('run Hidden Markov Model for (16S, 23S, 5S)')
		rrna.runHMM(fasta_file_ids, outLog=get_log_file_name(log_dir, 'hmm16S'))
		#print('run (16S, 23S, 5S) classification')
		#rrna.classify16S(fasta_file_ids, outLog=get_log_file_name(log_dir, 'mothurClassify16S'))
		#rrna.classify23S(fasta_file_ids, outLog=get_log_file_name(log_dir, 'mothurClassify23S'))
		#rrna.classify5S(fasta_file_ids, outLog=get_log_file_name(log_dir, 'mothurClassify5S'))

	# marker gene analysis (Amphora)
	if args.g:
		mg = MarkerGeneAnalysis(config, mg_database, working_dir, mg_working_dir)
		print 'run 31 marker gene analysis'
		mg.runMarkerGeneAnalysis(fasta_file_ids, outLog=get_log_file_name(log_dir, 'amphoraMG'))


def get_log_file_name(log_dir, description):
	now = datetime.datetime.now()
	time_stamp = str(str(now.year) + str(now.month) + str(now.day) + '_' + str(now.hour) + str(now.minute) + str(now.second))
	file_name = str(description) + '_' + time_stamp + '.txt'
	return os.path.normpath(os.path.join(log_dir, file_name))


if __name__ == "__main__":
	main()
