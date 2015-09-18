#!/usr/bin/python

__original_author__ = 'majda'
__author__ = 'peter hofmann'
__version__ = '0.0.3'


import sys
import os
import random
import argparse
import tempfile
import StringIO
from scripts.parallel import TaskCmd, runCmdParallel, reportFailedCmd
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.GenomePreparation.genomepreparation import GenomePreparation
# import maf_converter


class ReadSimulationWrapper(GenomePreparation):
	"""
	Default Class for all read simulation wrappers

	# TODO: validate genome: description still a problem for art illumina?
	"""
	_label = "ReadSimulationWrapper"

	def __init__(
		self, file_path_executable,
		separator='\t', max_processes=1, logfile=None, verbose=True, debug=False, seed=None, tmp_dir=None):
		"""
		Constructor

		@param file_path_executable:
		@type file_path_executable: str | unicode
		@param separator: separator to be expected in metadata files
		@type separator: str | unicode
		@param max_processes: Maximum number of processors simulating reads at the same time
		@type max_processes: int | long
		@param logfile: file handler or file path to a log file
		@type logfile: file | FileIO | StringIO | basestring
		@param verbose: Not verbose means that only warnings and errors will be past to stream
		@type verbose: bool
		@param debug: If true temporary files will be kept
		@type debug: bool
		@param seed: Seed used for read simulator, if option available
		@type seed: object
		@param tmp_dir: Directory for storage of temporary files
		@type tmp_dir: int | long
		"""
		assert self.validate_file(file_path_executable, executable=True)
		assert isinstance(separator, basestring)
		assert isinstance(max_processes, (int, long))
		assert isinstance(verbose, bool)
		assert isinstance(debug, bool)
		assert seed is None or isinstance(seed, (long, int, float, basestring))
		assert tmp_dir is None or isinstance(tmp_dir, basestring)
		if tmp_dir is not None:
			assert self.validate_dir(tmp_dir)
		else:
			tmp_dir = tempfile.gettempdir()
		self._tmp_dir = self.get_full_path(tmp_dir)
		self._debug = debug
		if seed is not None:
			random.seed(seed)
			# seed = abs(hash(seed))
			# assert len(str(seed)) > 4, "Seed '{}' is too short!".format(seed)
		# self._seed = abs(hash(seed))
		super(ReadSimulationWrapper, self).__init__(logfile=logfile, verbose=verbose)
		self._max_processes = max_processes
		self._separator = separator
		self._file_path_executable = file_path_executable
		self._fragments_size_mean = float(270)
		self._fragment_size_standard_deviation = int(round(self._fragments_size_mean * 0.1))
		self._read_length = 150
		self._temporary_files = set()

	def _close(self):
		"""
		Remove temporary files
		"""
		self._logger = None
		# delete temporary files
		self._remove_temporary_files()

	@staticmethod
	def _get_seed():
		return random.randint(0, sys.maxsize)

	def _remove_temporary_files(self):
		if self._debug:
			return
		while len(self._temporary_files) > 0:
			file_path = self._temporary_files.pop()
			if os.path.isfile(file_path):
				os.remove(file_path)

	# read genome location file
	def _read_genome_location_file(self, file_path):
		"""
		Read file with the file paths of gnomes

		@param file_path: File genome id associated with the file path of a genome
		@type file_path: str | unicode

		@return: Dictionary of genome id to file path
		@rtype: dict[str|unicode, str|unicode]
		"""
		self._logger.info('Reading genome location file')
		assert self.validate_file(file_path)
		dict_id_file_path = {}
		metadata_table = MetadataTable(logfile=self._logfile, verbose=self._verbose, separator=self._separator)
		iterator_distributions = metadata_table.parse_file(file_path, as_list=True)
		for genome_id, file_path_genome in iterator_distributions:
			assert genome_id != '', "Invalid genomid: '{}'".format(genome_id)
			assert file_path_genome != '', "Invalid file path: '{}'".format(genome_id)
			assert self.validate_file(file_path_genome), "Invalid file path: '{}'".format(genome_id)

			# check uniqueness
			assert genome_id not in dict_id_file_path, "Genome '{}' not unique in the distribution file!".format(genome_id)
			dict_id_file_path[genome_id] = file_path_genome
		return dict_id_file_path

	def _read_distribution_file(self, file_path):
		"""
		Read file with the distribution of a sample

		@param file_path: File genome id associated with the abundance of a genome
		@type file_path: str | unicode

		@return: Dictionary of genome id to file path
		@rtype: dict[str|unicode, float]
		"""
		self._logger.info('Reading distribution file')
		assert self.validate_file(file_path)
		dict_id_abundance = {}
		# dict_id_file_path = {}
		metadata_table = MetadataTable(logfile=self._logfile, verbose=self._verbose, separator=self._separator)
		iterator_distributions = metadata_table.parse_file(file_path, as_list=True)
		# for genome_id, abundance, genome_length, file_path_genome in iterator_distributions:
		for genome_id, abundance in iterator_distributions:
			assert genome_id != '', "Invalid genom id: '{}'".format(genome_id)
			assert abundance != '', "Invalid abundance: '{}'".format(genome_id)
			abundance = float(abundance)
			# assert file_path_genome != '', "Invalid file path: '{}'".format(genome_id)
			# assert self.validate_file(file_path_genome), "Invalid file path: '{}'".format(genome_id)
			assert self.validate_number(abundance, zero=False), "Invalid abundance: '{}'".format(genome_id)

			assert genome_id not in dict_id_abundance, "Genome '{}' not unique in the distribution file!".format(genome_id)
			# dict_id_file_path[genome_id] = file_path_genome
			dict_id_abundance[genome_id] = abundance
		# return dict_id_abundance, dict_id_file_path, relative_size_total, max_abundance
		return dict_id_abundance

	def get_multiplication_factor(
		self, dict_id_file_path, dict_id_abundance, total_size, min_sequence_length,
		file_format="fasta", sequence_type="dna", ambiguous=True):
		"""
		A factor is calculated based on total size of a sample to calculate the required covered
		# coverage = abundance * factor
		Files will be validated while the length of sequences are determined.

		@attention min_sequence_length: Sequences that are shorter than the expected fragment size are removed

		@param dict_id_file_path: Dictionary of genome id to file path
		@type dict_id_file_path: dict[str|unicode, str|unicode]
		@param dict_id_abundance: Dictionary of genome id to abundance
		@type dict_id_abundance: dict[str|unicode, float]
		@param total_size: Size of sample in base pairs
		@type total_size: int | long
		@param min_sequence_length: Minimum length of a sequence in base pairs
		@type min_sequence_length: int | long
		@param file_format: fasta or fastq
		@type file_format: str | unicode
		@param sequence_type: dna or rna or protein
		@type sequence_type: str | unicode
		@param ambiguous: DNA example for strict 'GATC',  ambiguous example 'GATCRYWSMKHBVDN'
		@type ambiguous: bool

		@return: Factor abundances will be multiplied by
		@rtype: float
		"""
		assert isinstance(dict_id_file_path, dict), "Expected dictionary, genome id as key, file path as value"
		assert isinstance(dict_id_abundance, dict), "Expected dictionary, genome id as key, abundance as value"
		assert isinstance(total_size, (int, long)), "Expected natural digit"
		assert isinstance(min_sequence_length, (int, long)), "Expected natural digit"
		assert isinstance(file_format, basestring), "Expected file format 'fasta'"
		assert isinstance(sequence_type, basestring), "Expected sequence type 'rna' or 'dna' or 'protein'"
		assert isinstance(ambiguous, bool)

		relative_size_total = 0
		for genome_id, abundance in dict_id_abundance.iteritems():
			try:
				min_seq_length, genome_length = self.get_sequence_lengths(
					file_path=dict_id_file_path[genome_id],
					file_format=file_format,
					sequence_type=sequence_type,
					ambiguous=ambiguous,
					key=None,
					silent=False)

				if min_seq_length < min_sequence_length:
					self._logger.info("Genome '{}' has sequences below minimum, creating filtered copy.".format(genome_id))
					new_file_path = self._remove_short_sequences(
						dict_id_file_path[genome_id], min_sequence_length, file_format="fasta")
					dict_id_file_path[genome_id] = new_file_path
					self._temporary_files.add(new_file_path)

			except IOError as e:
				self._remove_temporary_files()
				raise e

			relative_size = abundance * genome_length
			relative_size_total += relative_size
		return total_size / float(relative_size_total)

	def _remove_short_sequences(self, file_path, min_sequence_length, file_format="fasta"):
		"""
		Copies a genome with sequences shorter than a minimum removed.

		@param file_path: File genome id associated with the abundance of a genome
		@type file_path: str | unicode
		@param min_sequence_length: Minimum length of a sequence in base pairs
		@type min_sequence_length: int | long
		@param file_format:
		@type file_format: str | unicode

		@return: File path of the genome with removed short sequences
		@rtype: str | unicode
		"""
		assert self.validate_file(file_path)
		assert isinstance(min_sequence_length, (int, long)), "Expected natural digit"
		assert isinstance(file_format, basestring), "Expected file format 'fasta'"

		file_path_output = tempfile.mktemp(dir=self._tmp_dir)
		with open(file_path) as stream_input, open(file_path_output, 'w') as stream_output:
			total_base_pairs = self._stream_sequences_of_min_length(
				stream_input, stream_output,
				sequence_min_length=min_sequence_length,
				file_format=file_format
				)
			if total_base_pairs == 0:
				msg = "No valid sequences > {} found!".format(min_sequence_length)
				self._logger.error(msg)
				raise IOError(msg)
		return file_path_output


# #################
# ReadSimulationArt - Art-Illumina Wrapper
# #################


class ReadSimulationArt(ReadSimulationWrapper):
	"""
	Simulate reads using art illumina
	Currently pair-end reads only!
	"""
	_label = "ReadSimulationArtIllumina"

	_art_error_profiles = {
		"mi": "EmpMiSeq250R",
		"hi": "EmpHiSeq2kR",
		"hi150": "HiSeq2500L150R"}

	_art_read_length = {
		"mi": 250,
		"hi": 100,
		"hi150": 150}

	def __init__(self, file_path_executable, directory_error_profiles, **kwargs):
		super(ReadSimulationArt, self).__init__(file_path_executable, **kwargs)
		# check availability of profiles
		file_names_of_error_profiles = [
			filename+file_end
			for ep, filename in self._art_error_profiles.iteritems()
			for file_end in ['1.txt', '2.txt']
			]
		assert self.validate_dir(directory_error_profiles, file_names=file_names_of_error_profiles)
		# set default profile
		self._profile = "hi150"
		self._read_length = self._art_read_length["hi150"]
		self._directory_error_profiles = directory_error_profiles

	def simulate(
		self, file_path_distribution, file_path_genome_locations, directory_output,
		total_size, profile, fragments_size_mean, fragment_size_standard_deviation):
		"""
		Simulate reads based on a given sample distribution

		@param file_path_distribution: File genome id associated with the abundance of a genome
		@type file_path_distribution: str | unicode
		@param file_path_genome_locations: File genome id associated with the file path of a genome
		@type file_path_genome_locations: str | unicode
		@param directory_output: Directory for the sam and fastq files output
		@type directory_output: str | unicode
		@param total_size: Size of sample in base pairs
		@type total_size: int | long
		@param profile: Art illumina error profile: 'low', 'mi', 'hi', 'hi150'
		@type profile: str | unicode
		@param fragments_size_mean: Size of the fragment of which the ends are used as reads in base pairs
		@type fragments_size_mean: int | long
		@param fragment_size_standard_deviation: Standard deviation of the fragment size in base pairs.
		@type fragment_size_standard_deviation: int | long
		"""
		assert isinstance(total_size, (int, long)), "Expected natural digit"
		assert isinstance(fragments_size_mean, (int, long)), "Expected natural digit"
		assert isinstance(fragment_size_standard_deviation, (int, long)), "Expected natural digit"
		assert total_size > 0, "Total size needs to be a positive number"
		assert fragments_size_mean > 0, "Mean fragments size needs to be a positive number"
		assert fragment_size_standard_deviation > 0, "Fragment size standard deviation needs to be a positive number"
		assert self.validate_dir(directory_output)
		if profile is not None:
			assert profile in self._art_error_profiles, "Unknown art illumina profile: '{}'".format(profile)
			assert profile in self._art_read_length,  "Unknown art illumina profile: '{}'".format(profile)
			self._profile = profile
		if fragments_size_mean and fragment_size_standard_deviation:
			assert self.validate_number(fragments_size_mean, minimum=1)
			assert self.validate_number(fragment_size_standard_deviation, minimum=0)
			self._fragments_size_mean = fragments_size_mean
			self._fragment_size_standard_deviation = fragment_size_standard_deviation
		else:
			if fragment_size_standard_deviation:
				assert fragments_size_mean is not None, "Both, mean and sd are requires."
			if fragments_size_mean:
				assert fragment_size_standard_deviation is not None, "Both, mean and standard deviation, are required."
		self._logger.info("Using '{}' error profile.".format(profile))

		dict_id_abundance = self._read_distribution_file(file_path_distribution)
		dict_id_file_path = self._read_genome_location_file(file_path_genome_locations)
		assert set(dict_id_file_path.keys()).issuperset(dict_id_abundance.keys()), "Some ids do not have a genome location"

		min_sequence_length = self._fragments_size_mean - self._fragment_size_standard_deviation
		factor = self.get_multiplication_factor(
			dict_id_file_path, dict_id_abundance, total_size, min_sequence_length,
			file_format="fasta", sequence_type="dna", ambiguous=True)

		self._logger.debug("Multiplication factor: {}".format(factor))
		self._simulate_reads(dict_id_abundance, dict_id_file_path, factor, directory_output)

	# start ART readsimulator
	def _simulate_reads(self, dict_id_abundance, dict_id_file_path, factor, directory_output):
		"""
		Parallel simulation of reads

		@param dict_id_abundance: Dictionary of genome id to abundance
		@type dict_id_abundance: dict[str|unicode, float]
		@param dict_id_file_path: Dictionary of genome id to file path
		@type dict_id_file_path: dict[str|unicode, str|unicode]
		@param factor: Factor abundances will be multiplied by
		@type factor: float | int | long
		@param directory_output: Directory for the sam and fastq files output
		@type directory_output: str | unicode
		"""
		self._logger.info("Simulating reads using art Illumina readsimulator...")
		assert isinstance(dict_id_file_path, dict), "Expected dictionary, genome id as key, file path as value"
		assert isinstance(dict_id_abundance, dict), "Expected dictionary, genome id as key, abundance as value"
		assert isinstance(factor, (int, long, float)), "Factor must be a digit"
		assert self.validate_dir(directory_output)

		# add commands to a list of tasks to run them in parallel instead of calling them sequentially
		tasks = []
		for genome_id in dict_id_abundance.keys():
			file_path_input = dict_id_file_path[genome_id]
			abundance = dict_id_abundance[genome_id]
			fold_coverage = abundance * factor
			file_path_output_prefix = os.path.join(directory_output, str(genome_id))
			self._logger.debug("{id}\t{fold_coverage}".format(id=genome_id, fold_coverage=fold_coverage))
			system_command = self._get_sys_cmd(
				file_path_input=file_path_input,
				fold_coverage=fold_coverage,
				file_path_output_prefix=file_path_output_prefix)
			self._logger.debug("SysCmd: '{}'".format(system_command))
			self._logger.info("Simulating reads from {}: '{}'".format(genome_id, file_path_input))
			tasks.append(TaskCmd(system_command))
		list_of_fails = runCmdParallel(tasks, maxProc=self._max_processes)

		if list_of_fails is not None:
			self._logger.error("{} commands returned errors!".format(len(list_of_fails)))
			reportFailedCmd(list_of_fails)
		self._logger.info("Simulating reads finished")

	def _get_sys_cmd(self, file_path_input, fold_coverage, file_path_output_prefix):
		"""
		Build system command to be run.

		@param file_path_input: Path to genome fasta file
		@type file_path_input: str | unicode
		@param fold_coverage: coverage of a genome
		@type fold_coverage: int | long | float
		@param file_path_output_prefix: Output prefix used by art illumina
		@type file_path_output_prefix: str | unicode

		@return: System command to run art illumina
		@rtype: str | unicode
		"""
		assert self.validate_file(file_path_input)
		assert isinstance(fold_coverage, (int, long, float))
		assert self.validate_dir(file_path_output_prefix, only_parent=True)

		# TODO: mask 'N' default: '-nf 1'
		read_length = self._art_read_length[self._profile]
		error_profile = os.path.join(self._directory_error_profiles, self._art_error_profiles[self._profile])
		arguments = [
			"-sam", "-na",
			"-i '{}'".format(file_path_input),
			"-l", str(read_length),
			"-m", str(self._fragments_size_mean),
			"-s", str(self._fragment_size_standard_deviation),
			"-f", str(fold_coverage),
			"-o '{}'".format(file_path_output_prefix),
			"-1 '{}'".format(error_profile+'1.txt'),
			"-2 '{}'".format(error_profile+'2.txt'),
			]

		if self._logfile:
			arguments.append(">> '{}'".format(self._logfile))
		# else:
		# 	arguments.append("> /dev/null")

		# art illumina only accepts integer as seed!
		arguments.append("-rs '{}'".format(self._get_seed))

		cmd = "{exe} {args}".format(exe=self._file_path_executable, args=" ".join(arguments))
		return cmd


# #################
# #################
#
# FUTURE WORK:
# 	ReadSimulationPirs
# 	ReadSimulationPBSIM
#
# #################
# #################

# #################
# ReadSimulationPirs
# #################


# class ReadSimulationPirs(ReadSimulationWrapper):
# 	_label = "ReadSimulationPirs"
#
# 	def simulate(
# 		self, file_path_distributions, file_path_genome_locations, directory_output,
# 		total_size, read_length, fragments_size_mean, fragment_size_standard_deviation):
# 		raise Exception("Not fully implemented yet")
# 		assert self.validate_number(read_length, minimum=1)
# 		assert self.validate_number(fragments_size_mean, minimum=1)
# 		assert self.validate_number(fragment_size_standard_deviation, minimum=0)
# 		self._read_length = read_length
# 		self._fragments_size_mean = fragments_size_mean
# 		self._fragment_size_standard_deviation = fragment_size_standard_deviation
#
# 		dict_id_abundance = self._read_distribution_file(file_path_distributions)
# 		dict_id_file_path = self._read_genome_location_file(file_path_genome_locations)
#
# 		# coverage = abundance * factor
# 		# factor is calculated based on total size of sample
# 		min_sequence_length = self._fragments_size_mean - self._fragment_size_standard_deviation
# 		factor = self.get_multiplication_factor(
# 			dict_id_file_path, dict_id_abundance, total_size, min_sequence_length,
# 			file_format="fasta", sequence_type="dna", ambiguous=True)
# 		self._logger.debug("Multiplication factor: {}".format(factor))
# 		self._simulate_reads(dict_id_abundance, dict_id_file_path, factor, directory_output)
#
# 	# start pIRS readsimulator
# 	def _simulate_reads(self, dict_id_abundance, dict_id_file_path, factor, directory_output):
# 		# tmp_directory = "./"  # just write files locally
# 		# executable = os.path.join(self._directory_read_simulator, "pIRS/SimMetagenome.py")
# 		executable = self._file_path_executable
# 		# directory_temp = 'pIRS/nobackup/temp_abundance_file.csv'
# 		temp_abundance_filename = "temp_abundance_file.csv"
# 		for source_data_id in dict_id_abundance.keys():
# 			file_path_input = dict_id_file_path[source_data_id]
# 			abundance = dict_id_abundance[source_data_id]
# 			new_abundance = float(abundance) * factor
# 			with open(temp_abundance_filename, 'w') as temp_abundance_file:
# 				temp_abundance_file.write(source_data_id+'\t'+str(new_abundance)+'\n')
# 			with open(os.path.join(self._directory_read_simulator, 'pIRS/sampleConfig.cfg'), 'r') as config:
# 				with open("new_sampleConfig.cfg", 'w') as new_config:
# 					for line in config:
# 						if line.startswith("referenceSeq="):
# 							line = "referenceSeq=" + file_path_input + '\n'
# 						elif line.startswith('frequenciesInfo='):
# 							line = "frequenciesInfo=" + temp_abundance_filename + '\n'
# 						new_config.write(line)
# 				os.system("{} -c new_sampleConfig.cfg".format(executable))
#
#
# # #################
# # ReadSimulationPBSIM
# # #################
#
#
# class ReadSimulationPBSIM(ReadSimulationWrapper):
# 	_label = "ReadSimulationPBSIM"
#
# 	def simulate(
# 		self, file_path_distributions, file_path_genome_locations, directory_output,
# 		total_size, read_length, fragments_size_mean, fragment_size_standard_deviation):
# 		raise Exception("Not fully implemented yet")
# 		assert self.validate_number(read_length, minimum=1)
# 		assert self.validate_number(fragments_size_mean, minimum=1)
# 		assert self.validate_number(fragment_size_standard_deviation, minimum=0)
# 		self._read_length = read_length
# 		self._fragments_size_mean = fragments_size_mean
# 		self._fragment_size_standard_deviation = fragment_size_standard_deviation
#
# 		dict_id_abundance = self._read_distribution_file(file_path_distributions)
# 		dict_id_file_path = self._read_genome_location_file(file_path_genome_locations)
#
# 		# coverage = abundance * factor
# 		# factor is calculated based on total size of sample
# 		min_sequence_length = self._fragments_size_mean - self._fragment_size_standard_deviation
# 		factor = self.get_multiplication_factor(
# 			dict_id_file_path, dict_id_abundance, total_size, min_sequence_length,
# 			file_format="fasta", sequence_type="dna", ambiguous=True)
# 		self._logger.debug("Multiplication factor: {}".format(factor))
# 		self._simulate_reads(dict_id_abundance, dict_id_file_path, factor, directory_output)
#
# 	# start PBSIM readsimulator
# 	def _simulate_reads(self, dict_id_abundance, dict_id_file_path, factor, directory_output):
# 		# tmp_directory = "./"  # just write files locally
# 		# os.chdir(tmp_directory)
# 		# executable = os.path.join(self._directory_read_simulator, "pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim")
# 		executable = self._file_path_executable
# 		for source_data_id in dict_id_abundance.keys():
# 			file_path_input = dict_id_file_path[source_data_id]
# 			abundance = dict_id_abundance[source_data_id]
# 			# file_path_input, abundance, tax_id_predict, genome_name = dict_id_abundance[source_data_id]
# 			new_abundance = float(abundance) * factor
# 			# if not os.path.exists(tmp_directory+seq_id):
# 			#    os.mkdir(tmp_directory+seq_id)
# 			prefix = os.path.join(directory_output, str(source_data_id))
# 			arguments = [
# 				"--data-type", "CLR",
# 				"--depth", str(new_abundance),
# 				"--prefix", prefix,
# 				"--model_qc", os.path.join(self._directory_read_simulator, "pbsim-1.0.3-Linux-amd64/data/model_qc_clr"),
# 				file_path_input]
# 			# log_file.write("{} {}".format(executable, " ".join(arguments)))
# 			os.system("{} {} 2> /dev/null".format(executable, " ".join(arguments)))
#
# 		# if self._logger:
# 		#    self._logger.error("[ReadSimulation] pbsim currently not active. Please ask developer for more information.")
# 		# log_file.write(reader_folder+'pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim --data-type CLR
# 		#  --depth '+str(new_abundance)+'
# 		#  --model_qc '+reader_folder+'pbsim-1.0.3-Linux-amd64/data/model_qc_clr  '+address)
# 		# os.system(reader_folder+'pbsim-1.0.3-Linux-amd64/Linux-amd64/bin/pbsim --data-type CLR
# 		#  --depth '+str(new_abundance)+'
# 		#  --model_qc '+reader_folder+'pbsim-1.0.3-Linux-amd64/data/model_qc_clr  '+address)
# 		# TODO: file_path_output = tempfile.mktemp(dir=self._tmp_dir)
# 		maf_converter.main(directory_output)


# #################
# MAIN
# #################


def main(args=None):
	parser = argparse.ArgumentParser(description="Readsimulator")

	parser.add_argument(
		"total_size",
		default=None,
		type=long,
		help="Output size in base pairs, example 1GB: '1*10**9'")
	parser.add_argument(
		"-i", "--input_distribution",
		default=None,
		type=str,
		help="list with genome ID, abundance and genome size")
	parser.add_argument(
		"-l", "--input_genome_loaction",
		default=None,
		type=str,
		help="list with genome ID, and file path")
	parser.add_argument(
		"-o", "--directory_output",
		default=None,
		type=str,
		help="directory where fastq files will be written to")
	parser.add_argument(
		"-tmp", "--directory_tmp",
		default=None,
		type=str,
		help="temporary storage of files")
	parser.add_argument(
		"-log", "--file_path_logfile",
		default=None,
		type=str,
		help="Logfile")
	parser.add_argument(
		"-ep", "--error_profile",
		default="hi150",
		type=str,
		choices=["mi", "hi", "hi150"],
		help="mi: MiSeq, hi: HiSeq 100, hi150: HiSeq 150")
	parser.add_argument(
		"-s", "--fragment_size_standard_deviation", default=27, type=int,
		help="the standard deviation of DNA/RNA fragment size for paired-end simulations. (from art_sim)")
	parser.add_argument(
		"-m", "--fragments_size_mean", default=270, type=int,
		help="the mean size of DNA/RNA fragments for paired-end simulations. (from art_sim)")
	parser.add_argument(
		"-p", "--processor_pool", default=1, type=int,
		help="number of available processors")

	if args is None:
		options = parser.parse_args()
	else:
		options = parser.parse_args(args)

	directory_script = os.path.dirname(__file__)
	file_path_executable = os.path.join(directory_script, "tools", "readsimulator", "art_illumina")
	directory_error_profiles = os.path.join(directory_script, "tools", "readsimulator", "profile")

	simulator = ReadSimulationArt(
		file_path_executable=file_path_executable,
		directory_error_profiles=directory_error_profiles,
		separator='\t',
		max_processes=options.processor_pool,
		logfile=options.file_path_logfile,
		verbose=True,
		debug=False,
		seed=None,
		tmp_dir=options.directory_tmp)

	simulator.simulate(
		file_path_distribution=options.input_distribution,
		file_path_genome_locations=options.input_genome_loaction,
		directory_output=options.directory_output,
		total_size=options.total_size,
		profile=options.error_profile,
		fragments_size_mean=options.fragments_size_mean,
		fragment_size_standard_deviation=options.fragment_size_standard_deviation)


if __name__ == "__main__":
	main()
