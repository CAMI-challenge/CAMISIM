__author__ = 'hofmann'
__version__ = '0.0.4'


import sys
import os
import io
import errno
import itertools
import argparse
from Bio import SeqIO
from scripts.Validator.sequencevalidator import SequenceValidator


class FastaStreamer(SequenceValidator):
	"""
	Read files or list of files in a specific way and stream them to stdout
	Sequences are separated by '\000', a null character

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
			- if you have pairs or single blocks (pairs here: forward and backward reads)

	output is the same file with anonymous ids and a mapping form old ids to new ones

	# TODO: have user choose new_line character '\0', '\n'
	"""

	_legal_formats = ["fastq", "fasta"]

	def __init__(self, logfile=None, verbose=True, debug=False):
		super(FastaStreamer, self).__init__(logfile, verbose, debug, label="FastaStreamer")

	def stream_directory(self, directory, out_stream=sys.stdout, file_format="fastq", extension="fq", paired=False):
		"""
			Stream sequences consecutively

			@attention:

			@param directory: A directory
			@type directory: str
			@param out_stream: A stream the output will be written to.
			@type out_stream: file | io.FileIO | StringIO.StringIO
			@param file_format: Fasta format of input and output. Either 'fasta' or 'fastq'.
			@type file_format: str
			@param extension: file extension to be filtered for
			@type extension: str
			@param paired: sequences are streamed interweaved from a pair of files if True, else consecutively
			@type paired: bool

			@return: None
			@rtype: None
		"""
		assert isinstance(directory, str)
		directory = FastaStreamer.get_full_path(directory)
		assert self.validate_dir(directory)
		assert self.is_stream(out_stream)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in self._legal_formats
		assert extension is None or isinstance(extension, str)

		list_of_file = self.get_files_in_directory(directory, extension=extension)
		if not paired:
			self.consecutive_stream(list_of_file, out_stream=out_stream, file_format=file_format)
		else:
			self.interweave_stream(list_of_file, out_stream=out_stream, file_format=file_format, extension=extension)

	def stream_file(self, file_path, out_stream=sys.stdout, file_format="fastq", paired=False):
		"""
			Stream sequences consecutively

			@attention:

			@param file_path: A file path
			@type file_path: str
			@param out_stream: A stream the output will be written to.
			@type out_stream: file | io.FileIO | StringIO.StringIO
			@param file_format: Fasta format of input and output. Either 'fasta' or 'fastq'.
			@type file_format: str
			@param paired: sequences are streamed as pair, else one by one
			@type paired: bool

			@return: None
			@rtype: None
		"""
		assert isinstance(file_path, str)
		file_path = FastaStreamer.get_full_path(file_path)
		assert self.validate_file(file_path)
		assert self.is_stream(out_stream)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in self._legal_formats

		self.consecutive_stream(file_path, out_stream=out_stream, file_format=file_format, paired=paired)

	def consecutive_stream(self, src, out_stream=sys.stdout, file_format="fasta", paired=False):
		"""
			Stream sequences consecutively

			@attention:

			@param src: A file path or list of file paths
			@type src: str | list[str]
			@param out_stream: A stream the output will be written to.
			@type out_stream: file | io.FileIO | StringIO.StringIO
			@param file_format: Fasta format of input and output. Either 'fasta' or 'fastq'.
			@type file_format: str

			@return: None
			@rtype: None
		"""
		assert isinstance(src, (str, list))
		assert self.is_stream(out_stream)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in self._legal_formats

		list_of_file_paths = None
		if isinstance(src, str):
			assert self.validate_file(src)
			list_of_file_paths = [src]
		elif isinstance(src, list):
			for file_path in src:
				assert self.validate_file(file_path)
			list_of_file_paths = src

		for src in list_of_file_paths:
			if not os.path.exists(src):
				msg = "File does not exist: '{}'".format(src)
				self._logger.error(msg)
				raise IOError(msg)

		sequence_count = 0
		for src in list_of_file_paths:
			for seq_record in SeqIO.parse(src, file_format):
				new_line_suffix = '\0'
				if paired and sequence_count % 2 == 0:
					new_line_suffix = ""
				try:
					seq_record.description = ''
					out_stream.write(seq_record.format(file_format)+new_line_suffix)
				except IOError as e:
					if e.errno == errno.EPIPE or e.errno == errno.EINVAL:
						self._logger.error("Broken Pipe!")
						# Stop loop on "Invalid pipe" or "Invalid argument".
						# No sense in continuing with broken pipe.
						break
					else:
						self._logger.error("Something bad happened!!")
						# Raise any other error.
						raise
				sequence_count += 1

	# def interweave(self, directory=os.getcwd() + '/fastq/'):
	def interweave_stream(self, src, out_stream=sys.stdout, file_format="fasta", extension="fq"):
		"""
			Stream sequences of paired fasta files interweaved

			@attention:

			@param src: A file path or list of file paths
			@type src: str | list[str]
			@param out_stream: A stream the output will be written to.
			@type out_stream: file | io.FileIO | StringIO.StringIO
			@param file_format: Fasta format of input and output. Either 'fasta' or 'fastq'.
			@type file_format: str
			@param extension: file extension to be filtered for
			@type extension: str

			@return: None
			@rtype: None
		"""
		assert isinstance(src, (str, list))
		assert self.is_stream(out_stream)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in self._legal_formats
		assert isinstance(extension, str)

		list_of_file_paths = None
		if isinstance(src, str):
			assert self.validate_file(src)
			list_of_file_paths = [src]
		elif isinstance(src, list):
			for file_path in src:
				assert self.validate_file(file_path)
			list_of_file_paths = src

		file_path_one_suffix = "1.{}".format(extension)
		file_path_second_suffix = "2.{}".format(extension)
		for file_path in list_of_file_paths:
			if not file_path.endswith(file_path_one_suffix):
				continue
			file_path_one = file_path
			file_path_second = file_path[:file_path.rfind(file_path_one_suffix)] + file_path_second_suffix

			for seq_record_f, seq_record_b in itertools.zip_longest(SeqIO.parse(file_path_one, file_format), SeqIO.parse(file_path_second, file_format)):
				if seq_record_f is None or seq_record_b is None:
					msg = "forward and backward file have an unequal amount of sequences:\n"
					msg += "forward: '{}'\nbackward: '{}'\n".format(file_path_one, file_path_second)
					self._logger.error(msg)
					raise Exception(msg)
				out_stream.write(seq_record_f.format(file_format))
				out_stream.write(seq_record_b.format(file_format))
				out_stream.write('\0')

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
		help="file path or directory",
		action='store',
		type=str,
		default=None)
	parser.add_argument(
		"-out",
		help="output file, actually same as input file but with anonymous ids",
		action='store',
		type=argparse.FileType('w'),
		default=sys.stdout)
	parser.add_argument(
		"-s",
		help="process single sequences instead of pairs (e.g. single reads instead of paired end)",
		action="store_true",
		default=False)
	parser.add_argument(
		"-format",
		help="format of the input file fasta or fastq",
		action='store',
		choices=["fasta", "fastq"],
		default="fasta")
	parser.add_argument(
		"-ext",
		help="filter files in input directory using file extension",
		action='store',
		type=str,
		default=None)
	options = parser.parse_args()

	input_item = options.input
	file_format_chosen = options.format

	log_file_stream = options.log
	output_file_stream = options.out

	if input_item is None:
		sys.exit("[FastaStreamer] Bad input")

	fastastreamer = FastaStreamer(logfile=log_file_stream)
	if os.path.isfile(input_item):
		fastastreamer.stream_file(
			input_item,
			out_stream=output_file_stream,
			file_format=file_format_chosen,
			paired=not options.s
			)
	elif os.path.isdir(input_item):
		extension_chosen = options.ext
		fastastreamer.stream_directory(
			input_item,
			out_stream=output_file_stream,
			file_format=file_format_chosen,
			extension=extension_chosen,
			paired=not options.s
			)
	else:
		sys.exit("[FastaStreamer] Bad input: '{}'".format(input_item))

	if log_file_stream is not None:
		log_file_stream.close()
	if output_file_stream is not sys.stdout:
		output_file_stream.close()
