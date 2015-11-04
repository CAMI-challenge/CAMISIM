#!/usr/bin/env python

__author__ = 'hofmann'

import sys
import os
import argparse
from Bio import SeqIO
from scripts.Validator.sequencevalidator import SequenceValidator


def my_main():
	"""Parsing of arguments"""
	# example:
	# 	epilog = '''
	# ./concat_fasta_on_fasta.py \
	# 	-i input_file \
	# 	-o output_file \
	# '''

	description = "reads sequences from a fasta and concatinate it to another, renaming sequences based on file name"
	parser = argparse.ArgumentParser(description=description)

	parser.add_argument(
		"-i", "--input_file", default=None, type=str,
		help="path to file containing marker genes in multi fasta format")
	parser.add_argument(
		"-id", "--prefix_unique_id", default=None, type=str,
		help="unique id prefix of a genome, used as sequence name")
	parser.add_argument(
		"-o", "--output_file", default=None, type=str,
		help="target path to output file, all will be concatinated to this file")
	parser.add_argument(
		"-c", "--min_length", default=0, type=int,
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
		# parser.print_help()
		# sys.exit(1)
	else:
		if not os.path.isabs(output_file):
			# output_file = os.path.realpath(__file__)+"/"+output_file
			output_file = os.getcwd() + "/" + output_file
	
	if not os.path.exists(input_file):
		print "Error -i: input file does not exist"
		sys.exit(1)

	# with open(output_file, "a") as output_file_handle, open(out_bin_file, "a") as out_bin_file_handle:
	with open(output_file, "a") as output_file_handle:
		sequence_merger = SequenceMerger(
			stream_output=output_file_handle,
			stream_output_bin=None,
			stream_map_uid_sid=None,
			verbose=True)
		sequence_merger.merge(input_file, min_length=options.min_length, prefix_unique_id=options.unique_id)


class SequenceMerger(SequenceValidator):

	_label = "SequenceMerger"

	def __init__(self, stream_output, stream_output_bin, stream_map_uid_sid, logfile=None, verbose=False, debug=False):
		"""
		Constructor

		@param stream_output: Output stream
		@type stream_output: file | FileIO | StringIO
		@param stream_output_bin: Output stream of rejected sequences
		@type stream_output_bin: file | FileIO | StringIO
		@param stream_map_uid_sid: Output stream
		@type stream_map_uid_sid: file | FileIO | StringIO

		@rtype: None
		"""
		super(SequenceMerger, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		self._stream_output = stream_output
		self._stream_output_bin = stream_output_bin
		self._stream_map_uid_sid = stream_map_uid_sid

	def merge(self, file_path_input, min_length, prefix_unique_id=None, original_id=None, tax_id=""):
		"""
		Append sequences of fasta file to a stream

		@param file_path_input: File path to a fasta file with sequences
		@type file_path_input: str | unicode
		@param min_length: Minimum length of sequences
		@type min_length: int | long
		@param prefix_unique_id:
		@type prefix_unique_id: str | unicode

		@rtype: None
		"""
		assert isinstance(min_length, (int, long))
		assert isinstance(tax_id, basestring)
		assert isinstance(file_path_input, basestring)
		unique_id_set = set()
		# print input_file
		if prefix_unique_id is None:
			basename = os.path.basename(file_path_input)
			prefix_unique_id = os.path.splitext(basename)[0]

		if not os.path.exists(file_path_input):
			# sys.stderr.write("WARNING: [merge] File not found: '{file}'\n".format(file=input_file))
			self._logger.warning("No marker genes found for: '{unique_id}'".format(unique_id=prefix_unique_id))
			return

		counter = 0
		counter_rejected = 0
		for seq_record in SeqIO.parse(file_path_input, "fasta"):
			seq_length = len(seq_record.seq)
			if seq_record.id in unique_id_set:
				self._logger.warning("Removed duplicate entry of {}: {}".format(prefix_unique_id, seq_record.id))
				continue
			unique_id_set.add(seq_record.id)
			unique_sequence_id = "{}_{}".format(prefix_unique_id, counter)
			if seq_length >= min_length:
				counter += 1
				self._stream_output.write(">{}\n".format(unique_sequence_id))
				self._stream_output.writelines(seq_record.seq + "\n")
				if original_id is None:
					original_id = seq_record.id
				if self._stream_map_uid_sid and original_id:
					self._stream_map_uid_sid.write("{}\t{}\t{}\n".format(unique_sequence_id, original_id, tax_id))
			else:
				counter_rejected += 1
				# sys.stderr.write("WARNING: [merge] sequence too small. Size: {size} ID: '{uid}', File: '{file}'\n".format(
				# 				size=str(seq_length),
				# 				uid=unique_id,
				# 				file=os.path.basename(input_file)))
				if self._stream_output_bin:
					self._stream_output_bin.write(">{}\n".format(unique_sequence_id))
					self._stream_output_bin.writelines(seq_record.seq + "\n")

		if counter == 0:
			self._logger.warning("No valid marker gene found for {}: {}".format(prefix_unique_id, file_path_input))
		if counter_rejected > 0:
			self._logger.warning("{} marker genes rejected from {}".format(counter_rejected, prefix_unique_id))
		if counter > 0:
			self._logger.info("{} marker genes accepted from {}".format(counter, prefix_unique_id))


if __name__ == "__main__":
	# help(start_process)
	my_main()
