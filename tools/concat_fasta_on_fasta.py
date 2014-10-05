#!/usr/bin/env python

__author__ = 'hofmann'

import sys
import os
import argparse
from Bio import SeqIO

def my_main():
	"""Parsing of arguments"""
	# example:
	epilog = '''
./concat_fasta_on_fasta.py \
	-i input_file \
	-o output_file \
'''

	description = "reads sequences from a fasta and concatinate it to another, renaming sequences based on file name"
	parser = argparse.ArgumentParser(description=description)

	parser.add_argument("-i", "--input_file", default=None, type=str,
						help="path to file containing marker genes in multi fasta format")
	parser.add_argument("-id", "--unique_id", default=None, type=str,
						help="unique id of a genome, used as sequence name")
	parser.add_argument("-o", "--output_file", default=None, type=str,
						help="target path to output file, all will be concatinated to this file")
	parser.add_argument("-c", "--cutoff", default=0, type=int,
						help="minimal number of basepairs for a 16S gene")
	args = parser.parse_args()

	input_file = args.input_file
	if input_file is None:
		print "Error -i: Please pass a file of marker genes in multi fasta format"
		parser.print_help()
		sys.exit(1)

	output_file = args.output_file
	if output_file is None:
		print "Error -o: Please pass a file location for the output"
		#parser.print_help()
		#sys.exit(1)
	else:
		if not os.path.isabs(output_file):
			#output_file = os.path.realpath(__file__)+"/"+output_file
			output_file = os.getcwd()+"/"+output_file
	
	if not os.path.exists(input_file):
		print "Error -i: input file does not exist"
		sys.exit(1)
	
	#print input_file
	unique_id = args.unique_id
	if unique_id is None:
		temp = os.path.basename(input_file).split(".")
		unique_id = str(temp[0]) + "." + str(temp[1])
	
	counter = 0
	for seq_record in SeqIO.parse(input_file, "fasta"):
		seq_length = len(seq_record.seq)
		if seq_length >= args.cutoff:
			counter += 1
			with open(output_file, "a") as file_handler:
				file_handler.write(">{}_{}\n".format(unique_id, seq_record.id))
				file_handler.writelines(seq_record.seq + "\n")
		else:
			print "Warning: marker gene size below cutoff", "Size:", str(seq_length), "ID:", unique_id, "File:", os.path.basename(input_file)
			with open(output_file + "_below_cutoff.fna", "a") as file_handler:
				file_handler.write(">{}_{}\n".format(unique_id, seq_record.id))
				file_handler.writelines(seq_record.seq + "\n")

	if counter == 0:
		print "Warning: marker gene not found for {}: {}".format(unique_id, input_file)
	sys.exit(0)

if __name__ == "__main__":
	#help(start_process)
	my_main()


