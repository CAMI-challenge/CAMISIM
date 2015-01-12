#!/usr/bin/env python

__author__ = 'hofmann'

import sys
import os
import argparse
import subprocess

def start_process(input_file, output_file, config_file, processors):
	"""Starting the process of calculating the OTUs

	Required input: 16S fasta file, a alignments reference db and an threshold
	The fasta file should contain the 16S sequences of the unknown genomes and the reference genomes.
	The reference genomes will help finding a good threshold for spezies level, which should also distinguish known species from unknown

	input_file -- 16S multi fasta file
	alignments_reference_db -- GreenGenes/SILVA
	threshold -- Number between 0 and 1, default 0.03
	
	1. alignment of the 16S sequences
		mothur: 
			requires: 16S multi fasta file and alignments reference database file
			method: alignment ??improved by  refrence db??, filter emty columns, removing of not overlaping sequence tail and heads
			output: aligned, filtered multi fasta file
	
	2. calculation distance matrix
		options:
		a) phyml:
			requires: aligned multi fasta phylip formated
			method: BioNJ tree, with optimized branches, corrected dna distances
			output: newick tree
			- Rscript takes newicktree and builds dist. matrix (phylip format)
		b) raxmlHPC:
			requires: aligned multi fasta phylip formated
			method: maximum likelihood tree, corrected dna distances, best of #
			output: newick tree
			- Rscript takes newicktree and builds dist. matrix (phylip format)
		c) mothur:
			requires: aligned multi fasta
			method: uncorrected dna distance (like DNAdist)
			output: dist. matrix (phylip format)
		d) fasttree?

	3. clustering
		mothur
		requires: dist. matrix (phylip format) and distance threshold for clustering
		output: mothur formated OTU list
		"""
	subprocess.call(['./hmm', str(input_file), str(output_file), str(config_file), str(processors)])

def my_main():
	"""Parsing of arguments"""
	# example:
	epilog = '''
python run.py \
	-i input_file \
	-o output_file \
	-c config_file.cfg\
'''

	description = "Script to calculate otu from a marker gene alignment"
	parser = argparse.ArgumentParser(description=description)
	#parser = argparse.ArgumentParser(description=description, epilog=epilog)

	parser.add_argument("-i", "--input_file", default=None, type=str,
						help="path to file containing marker genes in multi fasta format")
	parser.add_argument("-o", "--output_file", default=None, type=str,
						help="target path to output file, will contain otu list in mothur format")
	parser.add_argument("-c", "--config_file", default=None, type=str,
						help="path to a config file")
	parser.add_argument("-p", "--processors", default=2, type=int,
						help="number of processors to be used")
	args = parser.parse_args()

	if args.input_file is None:
		print "Error -i: Please pass a file of marker genes in multi fasta format"
		parser.print_help()
		sys.exit(1)
	else:
		input_file = args.input_file

	if args.output_file is None:
		print "Error -o: Please pass a file location for the output"
		#parser.print_help()
		#sys.exit(1)
	else:
		output_file = args.output_file
		if not os.path.isabs(output_file):
			#output_file = os.path.realpath(__file__)+"/"+output_file
			output_file = os.getcwd()+"/"+output_file


	if args.config_file is None:
		print "Error -db: Please pass a file of markergenes in multi fasta format"
		parser.print_help()
		sys.exit(1)
	else:
		config_file = args.config_file

	if args.processors < 1:
		print "Error -p: only positiv number of processors legal"
		parser.print_help()
		sys.exit(1)
	else:
		processors = args.processors
	#print input_file
	start_process(input_file, output_file, config_file, processors)

	print "finished"
	sys.exit(0)


if __name__ == "__main__":
	#help(start_process)
	my_main()


