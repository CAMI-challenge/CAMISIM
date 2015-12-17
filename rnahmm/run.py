#!/usr/bin/env python

"""
	Copyright (C) 2014  Ivan Gregor

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

	Note that we could have written some parts of this code in a nicer way,
	but didn't have time. Be careful when reusing the source code.


	Master script of the PPSplus.
"""

import os
import argparse
from algbioi.com import common
from algbioi.com import taxonomy_ncbi
from algbioi.core.sequences import Sequences
from algbioi.core.analysis16s import RRNA16S
from scripts.configparserwrapper import ConfigParserWrapper


def main():
	# external commands will be executed in Shell in Unix/Linux
	assert os.name == 'posix', str('The pipeline runs only on "posix" systems (i.e. Unix/Linux compatible). ' +
		'Your system is "' + os.name + '"')

	parser = argparse.ArgumentParser(
		description='''PhyloPythiaS Plus is an extension of PhyloPythiaS.''',
		epilog='''Read the user documentation for more details.''')

	parser.add_argument(
		'-c', '--config',
		type=file, required=True,
		help='configuration file of the pipeline',
		metavar='config.cfg', dest='config')

	parser.add_argument(
		'-i', '--input_fasta_file', default=None, required=False, help='path to fasta file',
		metavar='fasta_file.fn', dest='i')

	parser.add_argument(
		'-out', '--output_folder', default=None, required=False, help='path to output folder',
		dest='out')

	parser.add_argument(
		'-n', '--runn-rrna16S', action='store_true',
		help='run hidden markov model searching for 16S, 23S, and 5S marker genes',
		dest='n')

	parser.add_argument(
		'-g', '--run-marker-gene-analysis', action='store_true',
		help='run hidden markov model and classify according to the "31" marker genes',
		dest='g')

	parser.add_argument(
		"-hmmer", "--hmmer", default=3, type=int,
		help="'2': rnammer; '3': hmmsearch using hmmer 3.0")

	parser.add_argument(
		"-log", "--logfile",
		type=str,
		default=None,
		help="pipeline output will written to this log file")

	args = parser.parse_args()

	# read configuration
	logfile = args.logfile
	config = ConfigParserWrapper(args.config, logfile=logfile, verbose=True)

	# pipeline directory
	pipeline_dir = args.out
	if pipeline_dir is None or not os.path.isdir(pipeline_dir):
		print("Pipeline directory doesn't exist: ", pipeline_dir)
		return

	# create the following directories in the working/output directory if they don't exist
	working_dir = os.path.join(pipeline_dir, 'working')
	mg_working_dir = os.path.join(working_dir, 'mgWorking')
	output_dir = os.path.join(pipeline_dir, 'output')
	dir_array = [
		working_dir, output_dir, os.path.join(working_dir, 'projectDir'),
		os.path.join(working_dir, 'sampleSpecificDir'), mg_working_dir]
	for dirPath in dir_array:
		if not os.path.isdir(dirPath):
			try:
				os.mkdir(dirPath)
			except OSError:
				print("Can't create directory", dirPath)
				return

	input_fasta_file = args.i
	if (input_fasta_file is None) or ((input_fasta_file is not None) and (not os.path.isfile(input_fasta_file))):
		print("The input fasta file %s doesn't exist" % input_fasta_file)
		return
	# read input fasta: contigs, scaffolds, mapping

	input_fasta_scaffolds_file = None
	scaffolds_to_contigs_map_file = None

	# create input id files (contigs, scaffolds, mappings)
	fasta_file_ids = common.createTagFilePath(working_dir, input_fasta_file, 'ids')
	seq_name_seq_id_file = common.createTagFilePath(working_dir, input_fasta_file, 'cToIds')

	scaffold_contig_map_ids_file = common.createTagFilePath(working_dir, input_fasta_file, 'mapSCIds')

	taxonomic_ranks = taxonomy_ncbi.TAXONOMIC_RANKS[1:]  # without root
	try:
		min_seq_len = config.get_value("MarkerGeneExtraction", "minSeqLen", is_digit=True)
	except Exception as e:
		print("Can't parse configuration entry (minSeqLen), make sure that it's an integer number")
		raise e

	# generates working fasta files always when the configuration file is newer than particular files to be generated
	sequences = Sequences(
		input_fasta_file, input_fasta_scaffolds_file, scaffolds_to_contigs_map_file, taxonomic_ranks, min_seq_len)
	if sequences.get_sequence_count() == -1:
		print("WARNING: [16S detector] File contains no valid sequences: {}".format(input_fasta_file))
		return

	sequences.writeSequences(fasta_file_ids)
	# print('Working contigs input fasta file created: %s' % fasta_file_ids)
	sequences.writeSeqNameSeqId(seq_name_seq_id_file)
	# print('Ids mapping for the working contigs fasta file created: %s' % seq_name_seq_id_file)
	sequences.writeScaffoldContigMap(scaffold_contig_map_ids_file)
	# print('Scaffolds -> contigs map ids file created: %s' % scaffold_contig_map_ids_file)
	# assert Common.seqFileCmp(inputFastaFile, fastaFileIds), 'The fasta IDs file contains different sequences!'

	# is it specified what to do?
	if not (args.n or args.g):
		print('Choose what do you want to do!')
		print(parser.print_help())
		return

	# run 16S analysis
	if args.n:
		hmmer = args.hmmer
		rrna = RRNA16S(config, None, working_dir)
		print('run Hidden Markov Model for (16S, 23S, 5S)')
		rrna.runHMM(fasta_file_ids, outLog=logfile, hmmer=hmmer, moltypes="ssu", kingdoms="arc,bac")


if __name__ == "__main__":
	main()
