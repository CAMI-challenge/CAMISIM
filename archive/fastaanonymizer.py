__author__ = 'cami'
__version__ = '0.0.8'

import sys
import os
import io
import random
import tempfile
import subprocess
from scripts.Validator.sequencevalidator import SequenceValidator


class FastaAnonymizer(SequenceValidator):
	"""
	Anonymization pipeline
	It takes a input file or directory, shuffles the sequences and replaces the sequence ids.
	It also can handle sequences that come in pairs

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

	_legal_formats = ["fastq", "fasta"]
	_random_source_file_path = None

	def __init__(self, logfile=None, verbose=True, debug=False, seed=None, tmp_dir=None):
		"""
			Anonymize fasta sequences

			@attention: 'shuf' is used which loads everything into memory!

			@param logfile: file handler or file path to a log file
			@type logfile: file | io.FileIO | StringIO.StringIO | str | unicode
			@param verbose: Not verbose means that only warnings and errors will be past to stream
			@type verbose: bool
			@param debug: more output and files are kept, manual clean up required
			@type debug: bool
			@param seed: The seed written to the random_source file used by the 'shuf' command
			@type seed: long | int | float | str | unicode
			@param tmp_dir: directory for temporary files, like the random_source file for 'shuf'
			@type tmp_dir: str | unicode

			@return: None
			@rtype: None
		"""
		assert isinstance(verbose, bool)
		assert isinstance(debug, bool)
		assert seed is None or isinstance(seed, (long, int, float, str))
		assert tmp_dir is None or isinstance(tmp_dir, str)
		if tmp_dir is not None:
			assert self.validate_dir(tmp_dir)
		else:
			tmp_dir = tempfile.gettempdir()
		self._tmp_dir = tmp_dir
		super(FastaAnonymizer, self).__init__(logfile, verbose, debug, label="FastaAnonymizer")

		if seed is not None:
			random.seed(seed)

		script_dir = os.path.dirname(self.get_full_path(__file__))
		self._anonymizer = os.path.join(script_dir, "anonymizer.py")
		self._fastastreamer = os.path.join(script_dir, "fastastreamer.py")
		assert self.validate_file(self._anonymizer)
		assert self.validate_file(self._fastastreamer)

	def _close(self):
		self._logger = None
		if self._debug:
			return

	@staticmethod
	def _get_seed():
		return random.randint(0, sys.maxsize)

	def get_command(
		self, file_path_mapping, path_input, file_path_output,
		sequence_prefix, file_format, paired=False, file_extension=None):
		"""
			Get the normalized absolute path.

			@attention:

			@param file_path_mapping: file path where the mapping should be saved
			@type file_path_mapping: str | unicode
			@param path_input: file path or directory of the input file(s)
			@type path_input: str | unicode
			@param file_path_output: output file path where the anonymous sequences should be saved
			@type file_path_output: str | unicode
			@param sequence_prefix: Prefix of the anonymous sequence id.
			@type sequence_prefix: str | unicode
			@param file_format: Fasta format of input and output. Either 'fasta' or 'fastq'.
			@type file_format: str | unicode
			@param paired: sequences are streamed as pair, else one by one
			@type paired: bool
			@param file_extension: file extension to be filtered for
			@type file_extension: str | unicode | None

			@return: System command line
			@rtype: str
		"""
		assert isinstance(path_input, str)
		assert self.validate_dir(path_input, silent=True) or self.validate_file(path_input, silent=True)
		assert isinstance(file_path_output, str)
		assert self.validate_dir(file_path_output, only_parent=True)
		assert isinstance(file_path_mapping, str)
		assert self.validate_dir(file_path_mapping, only_parent=True)
		assert isinstance(sequence_prefix, str)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in self._legal_formats
		assert file_extension is None or isinstance(file_extension, str)
		assert isinstance(paired, bool)

		# https://www.gnu.org/software/coreutils/manual/html_node/Random-sources.html#Random-sources
		continuous_random_byte_stream = \
			'get_seeded_random() { seed="$1"; openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt < /dev/zero 2>/dev/null; }; '

		# ##################
		# fastastreamer
		# ##################
		fastastreamer_args = [
			"-input '{}'".format(path_input),
			# "-out '{}'".format(file_path_output),
			"-format '{}'".format(file_format),
			]
		if file_extension is not None:
			fastastreamer_args.append("-ext '{}'".format(file_extension))

		# ##################
		# shuffle
		# ##################
		shuffle = ["shuf"]

		shuffle_args = ["-z"]
		seed = self._get_seed()
		shuffle_args.append("--random-source=<(get_seeded_random {seed})".format(seed=seed))

		# ##################
		# anonymizer args
		# ##################
		anonymizer_args = [
			"-prefix '{}'".format(sequence_prefix),
			"-format '{}'".format(file_format),
			"-map '{}'".format(file_path_mapping),
			"-out '{}'".format(file_path_output)]

		if not paired:
			anonymizer_args.append("-s")
			fastastreamer_args.append("-s")

		if self._logfile is not None:
			anonymizer_args.append("-log '{}'".format(self._logfile))
			fastastreamer_args.append("-log '{}'".format(self._logfile))

		command = continuous_random_byte_stream + "python3 '{fastastreamer}' {fastastreamer_args}".format(
			fastastreamer=self._fastastreamer,
			fastastreamer_args=" ".join(fastastreamer_args),
			)
		command += " | {shuf} {shuffle_args}".format(
			shuf=" ".join(shuffle),
			shuffle_args=" ".join(shuffle_args),
			)
		command += " | tr -d '\\000' | python3 '{anonymizer}' {anonymizer_args}".format(
			anonymizer=self._anonymizer,
			anonymizer_args=" ".join(anonymizer_args)
			)
		self._logger.debug(command)
		return command

	def shuffle_anonymize(
		self, path_input, file_path_output=None, file_path_mapping=None, prefix="", file_format=None, file_extension=None):
		"""
			Shuffle and anonymize sequences

			@attention:

			@param path_input: A directory or file path
			@type path_input: str | unicode
			@param file_path_output: Path of file the output will be written to.
			@type file_path_output: str | unicode
			@param file_path_mapping: Path of file the sequence id mapping will be written to.
			@type file_path_mapping: str | unicode
			@param prefix: Prefix of the anonymous sequence id.
			@type prefix: str | unicode
			@param file_format: Fasta format of input and output. Either 'fasta' or 'fastq'.
			@type file_format: str | unicode
			@param file_extension: file extension to be filtered for
			@type file_extension: str | unicode | None

			@return: file path of original to anonymous id mapping
			@rtype: str | unicode, str | unicode
		"""
		if file_path_output is None:
			file_path_output = tempfile.mktemp(dir=self._tmp_dir)
		if file_path_mapping is None:
			file_path_mapping = tempfile.mktemp(dir=self._tmp_dir)
		assert isinstance(path_input, str)
		assert self.validate_dir(path_input, silent=True) or self.validate_file(path_input, silent=True)
		assert isinstance(file_path_output, str)
		assert self.validate_dir(file_path_output, only_parent=True)
		assert isinstance(file_path_mapping, str)
		assert self.validate_dir(file_path_mapping, only_parent=True)
		assert isinstance(prefix, str)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in self._legal_formats

		self._logger.info("Shuffle and anonymize '{}'".format(path_input))

		command = self.get_command(
			file_path_mapping=file_path_mapping,
			path_input=path_input,
			file_path_output=file_path_output,
			sequence_prefix=prefix,
			file_format=file_format,
			paired=False,
			file_extension=file_extension)

		exit_status = subprocess.call(command, shell=True, executable="bash")
		if not exit_status == 0:
			msg = "Error occurred anonymizing '{}'".format(path_input)
			self._logger.error(msg)
			raise OSError(msg)
		return file_path_output, file_path_mapping

	def interweave_shuffle_anonymize(
		self, path_input, file_path_output=None, file_path_mapping=None, prefix="", file_format=None, file_extension=None):
		"""
			Interweave, shuffle and anonymize sequences

			@attention:

			@param path_input: A directory or file path
			@type path_input: str | unicode
			@param file_path_output: Path of file the output will be written to.
			@type file_path_output: str | unicode
			@param file_path_mapping: Path of file the sequence id mapping will be written to.
			@type file_path_mapping: str | unicode
			@param prefix: Prefix of the anonymous sequence id.
			@type prefix: str | unicode
			@param file_format: Fasta format of input and output. Either 'fasta' or 'fastq'.
			@type file_format: str | unicode
			@param file_extension: file extension to be filtered for
			@type file_extension: str | unicode | None

			@return: file path of original to anonymous id mapping
			@rtype: str | unicode, str | unicode
		"""
		if file_path_output is None:
			file_path_output = tempfile.mktemp(dir=self._tmp_dir)
		if file_path_mapping is None:
			file_path_mapping = tempfile.mktemp(dir=self._tmp_dir)
		assert isinstance(path_input, str)
		assert self.validate_dir(path_input, silent=True) or self.validate_file(path_input, silent=True)
		assert isinstance(file_path_output, str)
		assert self.validate_dir(file_path_output, only_parent=True)
		assert isinstance(file_path_mapping, str)
		assert self.validate_dir(file_path_mapping, only_parent=True)
		assert isinstance(prefix, str)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in self._legal_formats

		self._logger.info("Interweave shuffle and anonymize")

		command = self.get_command(
			file_path_mapping=file_path_mapping,
			path_input=path_input,
			file_path_output=file_path_output,
			sequence_prefix=prefix,
			file_format=file_format,
			paired=True,
			file_extension=file_extension)

		exit_status = subprocess.call(command, shell=True, executable="bash")
		if not exit_status == 0:
			msg = "Error occurred anonymizing '{}'".format(path_input)
			self._logger.error(msg)
			raise OSError(msg)
		return file_path_output, file_path_mapping
