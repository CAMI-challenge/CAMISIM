__author__ = 'cami'
__version__ = '0.0.5'

import sys
import io
import math
import argparse
from Bio import SeqIO
from scripts.Validator.sequencevalidator import SequenceValidator


class Anonymizer(SequenceValidator):
	"""
	input are streamed sequences which have identifier that need to be anonymous
	by default it reads from std.in

		for example a fastq file including paired end reads:
			@readid/1
			sequence
			+
			quality
			@readid/2
			sequence
			+
			quality
		you need to specify
			- a character that introduces a header line (@ here)
			- a block size (4 here)
			- if you have pairs or single blocks (pairs here: forward and backward reads)
			- optional a prefix (e.g. A|S1|R for reads from sample 1 in dataset A)

	output is the same file with anonymous ids and a mapping form old ids to new ones

	"""

	_label = "Anonymizer"
	_legal_formats = ["fastq", "fasta"]

	def anonymize_sequences(
		self, mapping, input_stream=sys.stdin, output_stream=sys.stdout,
		sequence_prefix='', file_format="fasta"):
		"""
			Anonymize sequences and write a mapping file

			@attention: only 'fasta' and 'fastq' supported

			@param mapping: output stream for sequence id mapping
			@type mapping: file | io.FileIO | StringIO.StringIO
			@param input_stream: Input stream of fasta format data
			@type input_stream: file | io.FileIO | StringIO.StringIO
			@param output_stream: Output stream of anonymous fasta format data
			@type output_stream: file | io.FileIO | StringIO.StringIO
			@param sequence_prefix: Prefix of the anonymous sequence id.
			@type sequence_prefix: str
			@param file_format: Fasta format of input and output. Either 'fasta' or 'fastq'.
			@type file_format: str

			@return: None
			@rtype: None
		"""
		assert self.is_stream(input_stream)
		assert self.is_stream(output_stream)
		assert self.is_stream(mapping)
		assert isinstance(sequence_prefix, str)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in self._legal_formats

		sequence_counter = 0
		for seq_record in SeqIO.parse(input_stream, file_format):
			anonymous_id = "{pre}{ct}".format(pre=sequence_prefix, ct=sequence_counter)
			map_line = "{oldid}\t{newid}\n".format(oldid=seq_record.id, newid=anonymous_id)
			mapping.write(map_line)
			seq_record.id = anonymous_id
			seq_record.description = ''
			output_stream.write(seq_record.format(file_format))
			sequence_counter += 1

	def anonymize_sequence_pairs(
		self, mapping, input_stream=sys.stdin, output_stream=sys.stdout,
		sequence_prefix='', file_format="fasta"):
		"""
			Anonymize pairs of sequences

			@attention: only 'fasta' and 'fastq' supported

			@param mapping: output stream for sequence id mapping
			@type mapping: file | io.FileIO | StringIO.StringIO | None
			@param input_stream: Input stream of fasta format data
			@type input_stream: file | io.FileIO | StringIO.StringIO | None
			@param output_stream: Output stream of anonymous fasta format data
			@type output_stream: file | io.FileIO | StringIO.StringIO | None
			@param sequence_prefix: Prefix of the anonymous sequence id.
			@type sequence_prefix: str
			@param file_format: Fasta format of input and output. Either 'fasta' or 'fastq'.
			@type file_format: str

			@return: None
			@rtype: None
		"""
		assert self.is_stream(input_stream)
		assert self.is_stream(output_stream)
		assert self.is_stream(mapping)
		assert isinstance(sequence_prefix, str)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in self._legal_formats

		sequence_counter = 0
		for seq_record in SeqIO.parse(input_stream, file_format):
			mod_value = sequence_counter % 2
			anonymous_id = "{pre}{ct}/{pair}".format(
				pre=sequence_prefix,
				ct=int(math.floor(sequence_counter/2.0)),
				pair=mod_value+1)
			map_line = "{oldid}\t{newid}\n".format(oldid=seq_record.id, newid=anonymous_id)
			mapping.write(map_line)
			seq_record.id = anonymous_id
			seq_record.description = ''
			output_stream.write(seq_record.format(file_format))
			sequence_counter += 1

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-log",
		help="log file",
		action='store',
		type=argparse.FileType('a'),
		default=None)
	parser.add_argument(
		"-input",
		help="input file (e.g. shuffled reads), reads from std.in by default",
		action='store',
		type=argparse.FileType('r'),
		default=sys.stdin)
	parser.add_argument(
		"-out",
		help="output file, actually same as input file but with anonymous ids",
		action='store',
		type=argparse.FileType('w'),
		default=sys.stdout)
	parser.add_argument(
		"-map",
		help="mapping file: old IDs -> anonymous IDs",
		action='store',
		type=argparse.FileType('w'),
		default=None)
	parser.add_argument(
		"-prefix",
		help="prefix for the anonymous ids",
		action='store',
		default="")
	parser.add_argument(
		"-s",
		help="process single blocks instead of pairs (e.g. single reads instead of paired end)",
		action="store_true",
		default=False)
	parser.add_argument(
		"-format",
		help="format of the input file fasta or fastq",
		action='store',
		choices=["fasta", "fastq"],
		default="fasta")
	options = parser.parse_args()

	input_file_stream = options.input

	map_file_stream = options.map
	output_file_stream = options.out
	log_file_stream = options.log

	anonymizer = Anonymizer(logfile=log_file_stream)
	if options.s:
		anonymizer.anonymize_sequences(
			mapping=map_file_stream,
			input_stream=input_file_stream,
			output_stream=output_file_stream,
			sequence_prefix=options.prefix,
			file_format=options.format)
	else:
		anonymizer.anonymize_sequence_pairs(
			mapping=map_file_stream,
			input_stream=input_file_stream,
			output_stream=output_file_stream,
			sequence_prefix=options.prefix,
			file_format=options.format)

	map_file_stream.close()
	if output_file_stream is not sys.stdout:
		output_file_stream.close()
	if input_file_stream is not sys.stdin:
		input_file_stream.close()
