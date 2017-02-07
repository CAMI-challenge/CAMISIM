__author__ = 'Peter Hofmann'
# original prototype:
#   http://armchairbiology.blogspot.de/2013/11/ani-are-you-okay-are-you-okay-ani.html
#   (c) The James Hutton Institute 2013
#   Author: Leighton Pritchard
#
#   Contact:
#   leighton.pritchard@hutton.ac.uk
#   GNU General Public License


import sys
import traceback
import tempfile
import shutil
import scripts.parallel as parallel
from scripts.Validator.validator import Validator
from scripts.MetaDataTable.metadatatable import MetadataTable


try:
	from Bio import SeqIO
except ImportError:
	SeqIO = None
	print "Biopython required for script, but not found (exiting)"
	sys.exit(1)


class ANIm(Validator):
	"""calculation average nucleotide identity"""
	def __init__(
		self, file_path_query_genomes_location, file_path_reference_genomes_location, file_path_reference_taxid_map,
		file_path_nucmer="nucmer", minimum_alignment=0.8, separator='\t', temp_directory=None, max_processors=1,
		logfile=None, verbose=False, debug=False):
		"""
		Constructor

		@param file_path_query_genomes_location:
		@type file_path_query_genomes_location: str|unicode
		@param file_path_reference_genomes_location:
		@type file_path_reference_genomes_location: str|unicode
		@param file_path_reference_taxid_map:
		@type file_path_reference_taxid_map: str|unicode
		@param file_path_nucmer:
		@type file_path_nucmer: str|unicode
		@param minimum_alignment:
		@type minimum_alignment: str|unicode|int|long|float
		@param separator:
		@type separator: str|unicode
		@param temp_directory:
		@type temp_directory: str|unicode
		@param max_processors:
		@type max_processors: int|long
		@param logfile: file handler or file path to a log file
		@type logfile: file | FileIO | StringIO | basestring
		@param verbose: Not verbose means that only warnings and errors will be past to stream
		@type verbose: bool
		@param debug: Display debug messages
		@type debug: bool

		@rtype: None
		"""
		assert self.validate_file(file_path_query_genomes_location)
		assert self.validate_file(file_path_reference_genomes_location)
		assert self.validate_file(file_path_reference_taxid_map)
		assert self.validate_file(file_path_nucmer, executable=True)
		assert temp_directory is None or self.validate_dir(temp_directory)
		assert isinstance(minimum_alignment, (int, float))
		assert self.validate_number(minimum_alignment, minimum=0, maximum=1)
		assert isinstance(separator, basestring)
		assert isinstance(max_processors, (int, long))
		assert self.validate_number(max_processors, minimum=1)
		super(ANIm, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		self._CUM_RETVALS = 0
		self._max_processors = max_processors
		self._file_path_nucmer = file_path_nucmer
		self._tmp_dir = temp_directory
		self._separator = separator
		if temp_directory is None:
			self._tmp_dir = tempfile.mkdtemp()
		else:
			self._tmp_dir = tempfile.mkdtemp(dir=temp_directory)
		self._cmd_lines = []
		data_table = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		data_table.read(file_path_query_genomes_location)
		self._query_gid_to_location = data_table.get_map(0, 1)
		data_table.read(file_path_reference_genomes_location)
		self._reference_gid_to_location = data_table.get_map(0, 1)
		data_table.read(file_path_reference_taxid_map)
		self._reference_gid_to_taxid = data_table.get_map(0, 1)
		self._total_lengths = {}
		self._minimum_alignment = minimum_alignment
		self._used_file_names = {}

	def __exit__(self, type, value, traceback):
		super(ANIm, self).__exit__(type, value, traceback)
		if not self._debug:
			shutil.rmtree(self._tmp_dir)

	# def get_total_organism_length(self):
	# 	self.get_organism_lengths(self._reference_gid_to_location)
	# 	self.get_organism_lengths(self._query_gid_to_location)

	# Construct a command-line for NUCmer reference_id, candidate_id
	def get_nucmer_cmd(self, reference_id, candidate_id, mum=True, maxmatch=False):
		"""
		Construct a command-line for NUCmer pairwise comparison, and return as
		a string

		- f1, f2 are the locations of two input FASTA files for analysis

		We use the -mum option so that we consider matches that are unique in
		both the reference and the query. -mumreference gives us matches
		unique only in the reference and -maxmatch gives us matches to all
		regions, regardless of uniqueness. We may want to make this an option.

		@param reference_id:
		@type reference_id: str|unicode
		@param candidate_id:
		@type candidate_id: str|unicode
		@param mum:
		@type mum: bool
		@param maxmatch:
		@type maxmatch: bool

		@return: Command line
		@rtype: str|unicode
		"""
		assert isinstance(reference_id, basestring)
		assert isinstance(candidate_id, basestring)
		assert isinstance(mum, bool)
		assert isinstance(maxmatch, bool)
		# if reference_id not in self._total_lengths:
		# 	self._total_lengths[reference_id] = self.get_organism_length(self._reference_gid_to_location[reference_id])
		if candidate_id not in self._total_lengths:
			self._total_lengths[candidate_id] = self.get_organism_length(self._query_gid_to_location[candidate_id])
		out_file_name = tempfile.mktemp(dir=self._tmp_dir, prefix=str(len(self._used_file_names)))
		self._used_file_names[(reference_id, candidate_id)] = out_file_name + ".delta"
		# out_file_name = "{}_vs_{}".format(queri_id, ref_id)
		# prefix = os.path.join(self._tmp_dir, out_file_name)
		# Do we use the --maxmatch option?
		mode = ""
		if maxmatch:
			mode = "-maxmatch"
		elif mum:
			mode = "-mum"

		cmd = "{nucmer} {mode} -p {prefix} {reference} {query}".format(
			nucmer=self._file_path_nucmer, mode=mode, prefix=out_file_name,
			reference=self._reference_gid_to_location[reference_id],
			query=self._query_gid_to_location[candidate_id])
		return cmd

	# Run NUCmer pairwise on the input files, using multiprocessing
	def add_nucmer_cmd_lines(self, reference_ids, candidate_ids):
		"""
		We loop over all FASTA files in reference_file_names, generating NUCmer
		command lines for each reference_file_names comparison

		@param reference_ids:
		@type reference_ids: list[str|unicode]
		@param candidate_ids:
		@type candidate_ids: list[str|unicode]

		@rtype: None
		"""
		assert isinstance(reference_ids, list)
		assert isinstance(candidate_ids, list)
		# self._logger.info("add_nucmer_cmd_lines: %s %s" % (candidate_ids, reference_ids))
		unique_set = set()
		# print "C", candidate_ids
		# print ""
		for candidate_id in candidate_ids:
			# print "R", reference_ids
			# print ""
			# print "RF", self.reference_file_names.keys()
			for reference_id in reference_ids:
				tupel = (candidate_id, reference_id)
				if tupel in unique_set:
					continue
				unique_set.add(tupel)
				if reference_id in self._reference_gid_to_location:
					# self.cmd_lines.extend([self.get_nucmer_cmd(reference_id, candidate_id) for reference_id in reference_ids])
					self._cmd_lines.append(self.get_nucmer_cmd(reference_id, candidate_id))
				# else:
				# 	self._logger.warning("No genom for reference: {}".format(reference_id))

	# Run a set of command lines using multiprocessing
	def multiprocessing_run(self):
		"""
		Distributes the passed command-line jobs using multiprocessing.

		@rtype: None
		"""
		self._logger.info("Running {} jobs with multiprocessing".format(len(self._cmd_lines)))
		list_cmd_task = [parallel.TaskCmd(cmd, self._tmp_dir) for cmd in self._cmd_lines]
		fail_list = parallel.runCmdParallel(list_cmd_task, self._max_processors)
		if fail_list is not None:
			parallel.reportFailedCmd(fail_list)
			self._CUM_RETVALS = -1 * len(fail_list)
		self._logger.info("Multiprocessing jobs completed")

	# Run NUCmer pairwise on the input files, using multiprocessing
	def run_cmd_lines(self):
		"""
		Run NUCmer to generate pairwise alignment data for each of the
		input FASTA files.

		- filenames is a list of input FASTA filenames, from which NUCmer
			command lines are constructed

		We loop over all FASTA files in the input directory, generating NUCmer
		command lines for each pairwise comparison, and then pass those
		command lines to be run using multiprocessing.

		@rtype: None
		"""
		if len(self._cmd_lines) < 1:
			self._logger.warning("NUCmer command lines: No lines!")
			self._logger.debug("Ref {} / {}".format(len(self._reference_gid_to_location), len(self._query_gid_to_location)))
			return
		# self._logger.info("NUCmer command lines:\n\t%s" % '\n\t'.join(self.cmd_lines[0]))
		self._logger.debug("NUCmer command lines:\n\t{}".format(self._cmd_lines[0]))
		# return
		self.multiprocessing_run()
		if 0 < self._CUM_RETVALS:
			self._logger.error("At least one NUCmer comparison failed. ANIm may fail.")

	# Get lengths of sequence for each organism
	@staticmethod
	def get_organism_length(file_path_fasta):
		"""
		Get the total length of a specific fasta file

		@param file_path_fasta:
		@type file_path_fasta: str|unicode

		@rtype: int|long
		"""
		assert isinstance(file_path_fasta, basestring)
		return sum([len(s) for s in SeqIO.parse(file_path_fasta, 'fasta')])

	# Get lengths of sequence for each organism
	# def get_organism_lengths(self, organism_file_names):
	# 	""" Returns a dictionary of total input sequence lengths, keyed by
	# 		organism.
	#
	# 		Biopython's SeqIO module is used to parse all sequences in the FASTA
	# 		file corresponding to each organism, and the total base count in each
	# 		is obtained.
	#
	# 		NOTE: ambiguity symbols are not discounted.
	# 	"""
	# 	self._logger.info("Processing organism sequence lengths")
	# 	for organism in organism_file_names:
	# 		self._total_lengths[organism] = sum([len(s) for s in SeqIO.parse(organism_file_names[organism], 'fasta')])

	# Parse NUCmer delta file to get total alignment length and total sim_errors
	@staticmethod
	def parse_delta(file_path):
		"""
		Reads a NUCmer output .delta file, extracting the aligned length and
		number of similarity errors for each aligned uniquely-matched region,
		and returns the cumulative total for each as a tuple.

		- filename is the path to the input .delta file

		@param file_path:
		@type file_path: str|unicode

		@return: aln_length, sim_errors
		@rtype: tuple[int|long, int|long]
		"""
		assert isinstance(file_path, basestring)
		aln_length, sim_errors = 0, 0
		for line in [l.strip().split() for l in open(file_path, 'rU').readlines()]:
			if line[0] == 'NUCMER' or line[0].startswith('>'):  # Skip headers
				continue
			# We only want lines with seven columns:
			if len(line) == 7:
				aln_length += abs(int(line[1]) - int(line[0]))
				sim_errors += int(line[4])
		return aln_length, sim_errors

	# Report last exception as string
	@staticmethod
	def last_exception():
		"""
		Returns last exception as a string, or use in logging.

		@rtype: str|unicode
		"""
		exc_type, exc_value, exc_traceback = sys.exc_info()
		return ''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))

	# Parse NUCmer delta output to store alignment total length, sim_error,
	# and percentage identity, for each pairwise comparison
	def process_delta(self):
		"""
		Returns a tuple containing a list and four dictionaries. The list
		describes the names of all organisms (derived from their filenames).
		The dictionaries describe results for pairwise comparisons: total
		aligned lengths; similarity errors in those alignments; the percentage
		of aligned length that matches (ANIm); and the percentage of the
		pairwise comparison that is aligned.

		For the total aligned length, similarity error, and ANIm dictionaries,
		as these are triangular/symmetrical matrices we only key them by
		(query, subject), but as the percentage aligned measure depends on the
		sequence we calculate it against, we report (query, subject) and
		(subject, query) values.

		- org_lengths is a dictionary of total sequence lengths for each
			input sequence

		@return: lengths, sim_errors, perc_ids, perc_aln
		@rtype: tuple[dict[tuple[], int|long], dict[tuple[], int|long], dict[tuple[], float], dict[tuple[], float]]
		"""
		# delta_files = self.get_input_files('.delta')
		self._logger.info("Processing .delta files")
		# We store pairwise comparison lengths in dictionaries, keyed by organism
		# ID pair tuples:
		# perc_aln is useful, as it is a matrix of the minimum percentage of an
		# organism's genome involved in a pairwise alignment
		lengths, sim_errors, perc_ids, perc_aln = {}, {}, {}, {}
		for (reference_id, candidate_id), delta_filename in self._used_file_names.iteritems():
			self._logger.info("Query organism: %s; Reference organism: %s" % (candidate_id, reference_id))
			tot_length, tot_sim_error = self.parse_delta(delta_filename)
			perc_id = 0
			try:
				perc_id = 1 - 1. * tot_sim_error / tot_length
			except ZeroDivisionError:
				# If this is thrown, the proximate cause is an empty NUCmer output
				# The root cause may be a failed NUCmer run (when CUM_RETVALS>0)
				# or that one or more of the sequences is too distant for NUCmer
				# to identify a similarity.
				self._logger.error("One or more of the NUCmer output files contains no useable output.")
				if 0 < self._CUM_RETVALS:
					self._logger.error("One or more NUCmer runs failed. Please investigate.")
					self._logger.error("Please retry the NUCmer comparison for %s vs %s manually" % (candidate_id, reference_id))
				else:
					self._logger.error(
						"The NUCmer comparison between %s and %s " % (candidate_id, reference_id) +
						"has no usable output. The comparison may be " +
						"too distant for use. Consider using --maxmatch.")
				self._logger.error(self.last_exception())
			lengths[(candidate_id, reference_id)] = tot_length
			sim_errors[(candidate_id, reference_id)] = tot_sim_error
			perc_ids[(candidate_id, reference_id)] = perc_id
			perc_aln[(candidate_id, reference_id)] = 1. * tot_length / self._total_lengths[candidate_id]
		return lengths, sim_errors, perc_ids, perc_aln

	# METHOD: ANIm
	# This method uses NUCmer to calculate pairwise alignments for the input
	# organisms, without chopping sequences into fragments. We follow the method
	# of Richter et al. (2009)
	def calculate_anim(self):
		"""
		Calculate ANI by the ANIm method, as described in Richter et al (2009)
		Proc Natl Acad Sci USA 106: 19126-19131 doi:10.1073/pnas.0906412106.

		All FASTA format files (selected by suffix) in the input directory
		are compared against each other, pairwise, using NUCmer (which must
		be in the path). NUCmer output is stored in the output directory.

		The NUCmer .delta file output is parsed to obtain an alignment length
		and similarity error count for every unique region alignment between
		the two organisms, as represented by the sequences in the FASTA files.

		These are processed to give matrices of aligned sequence lengths,
		similarity error counts, average nucleotide identity (ANI) percentages,
		and minimum aligned percentage (of whole genome) for each pairwise
		comparison.

		The matrices are written to file in a plain text tab-separated format.

		@return: lengths, sim_errors, perc_ids, perc_aln
		@rtype: tuple[dict[tuple[], int|long], dict[tuple[], int|long], dict[tuple[], float], dict[tuple[], float]]
		"""
		self._logger.info("Running ANIm method")
		self.run_cmd_lines()
		# lengths, sim_errors, perc_ids, perc_aln = self.process_delta()
		return self.process_delta()

	def calculate_best_anim(self):
		"""
		Return the results of the best hit in the references

		@return: lengths, sim_errors, perc_ids, perc_aln
		@rtype: dict[str|unicode, int|long], dict[str|unicode, int|long], dict[str|unicode, float], dict[str|unicode, float]
		"""
		min_lengths, min_sim_errors, min_perc_ids, min_perc_aln, ncbi, best_ani_by_candidate_id = {}, {}, {}, {}, {}, {}
		lengths, sim_errors, perc_ids, perc_aln = self.calculate_anim()

		for candidate_id_reference_id in perc_ids:
			candidate_id = candidate_id_reference_id[0]
			reference_id = candidate_id_reference_id[1]
			# self._logger.info("{}: candidate_id: {}; reference_id: {}".format(
			# candidate_id_reference_id, candidate_id, reference_id))
			# ani_ish = perc_ids[candidate_id_reference_id] * perc_aln[candidate_id_reference_id]
			if perc_aln[candidate_id_reference_id] < self._minimum_alignment:
				continue
			if (candidate_id not in best_ani_by_candidate_id or
						perc_ids[candidate_id_reference_id] > best_ani_by_candidate_id[candidate_id]):
				best_ani_by_candidate_id[candidate_id] = perc_ids[candidate_id_reference_id]
				min_lengths[candidate_id] = lengths[candidate_id_reference_id]
				min_sim_errors[candidate_id] = sim_errors[candidate_id_reference_id]
				min_perc_ids[candidate_id] = perc_ids[candidate_id_reference_id]
				min_perc_aln[candidate_id] = perc_aln[candidate_id_reference_id]
				ncbi[candidate_id] = self._reference_gid_to_taxid[reference_id]  # reference_id.split('.')[0]
		return min_lengths, min_sim_errors, min_perc_ids, min_perc_aln, ncbi
