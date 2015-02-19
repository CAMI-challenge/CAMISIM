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
	parser.add_argument("-c", "--min_length", default=0, type=int,
						help="minimal number of basepairs for a 16S gene")
	options = parser.parse_args()

	input_file = options.input_file
	if input_file is None:
		print "Error -i: Please pass a file of marker genes in multi fasta format"
		parser.print_help()
		sys.exit(1)

	output_file = options.output_file
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

	merge(input_file, output_file, options.min_length, unique_id=options.unique_id)


def merge(input_file, output_file, min_length, unique_id=None, out_bin_file=None):
	unique_id_set = set()
	#print input_file
	if unique_id is None:
		basename = os.path.basename(input_file)
		unique_id = os.path.splitext(basename)[0]

	if not os.path.exists(input_file):
		#sys.stderr.write("WARNING: [merge] File not found: '{file}'\n".format(file=input_file))
		sys.stderr.write("WARNING: [merge] No marker genes found for: '{unique_id}'\n".format(unique_id=unique_id))
		return

	counter = 0
	counter_small = 0
	for seq_record in SeqIO.parse(input_file, "fasta"):
		seq_length = len(seq_record.seq)
		if seq_record.id in unique_id_set:
			sys.stderr.write("WARNING: [merge] Removed duplicate entry of {}: {}\n".format(unique_id, seq_record.id))
			continue
		unique_id_set.add(seq_record.id)
		if seq_length >= min_length:
			counter += 1
			with open(output_file, "a") as file_handler:
				file_handler.write(">{}_{}\n".format(unique_id, seq_record.id))
				file_handler.writelines(seq_record.seq + "\n")
		else:
			counter_small += 1
			#sys.stderr.write("WARNING: [merge] sequence too small. Size: {size} ID: '{uid}', File: '{file}'\n".format(
			#				size=str(seq_length),
			#				uid=unique_id,
			#				file=os.path.basename(input_file)))
			if out_bin_file:
				with open(out_bin_file, "a") as file_handler:
					file_handler.write(">{}_{}\n".format(unique_id, seq_record.id))
					file_handler.writelines(seq_record.seq + "\n")

	if counter == 0:
		sys.stderr.write("WARNING: [merge] No valid marker gene found for {}: {}\n".format(unique_id, input_file))
	if counter_small > 0:
		sys.stderr.write("WARNING: [merge] {} marker genes rejected from {}\n".format(counter_small, unique_id))

if __name__ == "__main__":
	#help(start_process)
	my_main()


