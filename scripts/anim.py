__author__ = 'Peter Hofmann'
# original prototype:
#   http://armchairbiology.blogspot.de/2013/11/ani-are-you-okay-are-you-okay-ani.html
#   (c) The James Hutton Institute 2013
#   Author: Leighton Pritchard
#
#   Contact:
#   leighton.pritchard@hutton.ac.uk
#   GNU General Public License


import os
import sys
import traceback
import tempfile
import shutil
import parallel
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
		# self.used_file_names = {}

	def __exit__(self, type, value, traceback):
		super(ANIm, self).__exit__(type, value, traceback)
		if not self._debug:
			shutil.rmtree(self._tmp_dir)

	def get_total_organism_length(self):
		self.get_organism_lengths(self._reference_gid_to_location)
		self.get_organism_lengths(self._query_gid_to_location)

	# Construct a command-line for NUCmer reference_id, candidate_id
	def get_nucmer_cmd(self, reference_id, candidate_id, mum=True, maxmatch=False):
		""" Construct a command-line for NUCmer pairwise comparison, and return as
			a string

			- f1, f2 are the locations of two input FASTA files for analysis

			We use the -mum option so that we consider matches that are unique in
			both the reference and the query. -mumreference gives us matches
			unique only in the reference and -maxmatch gives us matches to all
			regions, regardless of uniqueness. We may want to make this an option.
		"""
		if reference_id not in self._total_lengths:
			self._total_lengths[reference_id] = self.get_organism_length(self._reference_gid_to_location[reference_id])
		if candidate_id not in self._total_lengths:
			self._total_lengths[candidate_id] = self.get_organism_length(self._query_gid_to_location[candidate_id])
		if self._total_lengths[reference_id] < self._total_lengths[candidate_id]:
			ref_id = reference_id
			queri_id = candidate_id
		else:
			ref_id = candidate_id
			queri_id = reference_id
		# TODO: not from filename!!
		out_file_name = "{}_vs_{}".format(queri_id, ref_id)
		prefix = os.path.join(self._tmp_dir, out_file_name)
		# Do we use the --maxmatch option?
		mode = ""
		if maxmatch:
			mode = "-maxmatch"
		elif mum:
			mode = "-mum"

		cmd = "{nucmer} {mode} -p {prefix} {reference} {query}".format(
			nucmer=self._file_path_nucmer, mode=mode, prefix=prefix,
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

	# Multiprocessing callback to logger
	def logger_callback(self, val):
		""" Basic callback for multiprocessing just to log status of each job

			- val is an integer returned by multiprocessing, describing the run
				status
		"""
		self._logger.info("Multiprocessing run completed with status: %s" % val)
		# Keep track of returned values, as these help diagnose problems
		# for ANIm analyses
		self._CUM_RETVALS += val

	# Run a set of command lines using multiprocessing
	def multiprocessing_run(self):
		""" Distributes the passed command-line jobs using multiprocessing.

			- cmdlines is an iterable of command line strings
		"""
		self._logger.info("Running {} jobs with multiprocessing".format(len(self._cmd_lines)))
		# return
		cmd_task_list = [parallel.TaskCmd(cmd, self._tmp_dir) for cmd in self._cmd_lines]
		fail_list = parallel.runCmdParallel(cmd_task_list, self._max_processors)
		if fail_list is not None:
			parallel.reportFailedCmd(fail_list)
			self._CUM_RETVALS = -1 * len(fail_list)
		# pool = multiprocessing.Pool(self.pool_size)
		# completed = []
		# if verbose:
		# 	callback_fn = self._logger_callback
		# else:
		# 	callback_fn = completed.append
		# [pool.apply_async(subprocess.call, (str(cmd_line), ), {'stderr': subprocess.PIPE, 'shell': sys.platform != "win32"}, callback=callback_fn)
		# 	for cmd_line in self.cmd_lines]
		# pool.close()        # Run jobs
		# pool.join()         # Collect output
		# self._logger.info("Multiprocessing jobs completed:\n%s" % completed)
		self._logger.info("Multiprocessing jobs completed")

	# Run NUCmer pairwise on the input files, using multiprocessing
	def run_cmd_lines(self):
		""" Run NUCmer to generate pairwise alignment data for each of the
			input FASTA files.

			- filenames is a list of input FASTA filenames, from which NUCmer
				command lines are constructed

			We loop over all FASTA files in the input directory, generating NUCmer
			command lines for each pairwise comparison, and then pass those
			command lines to be run using multiprocessing.
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
	def get_organism_length(organism_file_name):
		"""
		"""
		return sum([len(s) for s in SeqIO.parse(organism_file_name, 'fasta')])

	# Get lengths of sequence for each organism
	def get_organism_lengths(self, organism_file_names):
		""" Returns a dictionary of total input sequence lengths, keyed by
			organism.

			Biopython's SeqIO module is used to parse all sequences in the FASTA
			file corresponding to each organism, and the total base count in each
			is obtained.

			NOTE: ambiguity symbols are not discounted.
		"""
		self._logger.info("Processing organism sequence lengths")
		for organism in organism_file_names:
			self._total_lengths[organism] = sum([len(s) for s in SeqIO.parse(organism_file_names[organism], 'fasta')])

	# Get filenames for each organism
	def get_organism_file_names(self, file_filepath_list):
		""" Returns a dictionary of filenames, keyed by organism.

			Biopython's SeqIO module is used to parse all sequences in the FASTA
			file corresponding to each organism, and the total base count in each
			is obtained.

			NOTE: ambiguity symbols are not discounted.
		"""
		self._logger.info("Processing organism filenames")
		file_names = {}
		with open(file_filepath_list, 'r') as file_handler:
			for line in file_handler:
				if len(line) > 0 and line[0] == "#":
					continue
				data = line.strip().split("\t")
				organism = data[0]
				filename = data[1]
				file_names[organism] = filename
		return file_names

	# Get list of FASTA files in a directory
	def get_input_files(self, *ext):
		""" Returns a list of files in the input directory with the passed
			extension

			- dir is the location of the directory containing the input files

			- *ext is a list of arguments describing permissible file extensions
		"""
		file_list = [filename for filename in os.listdir(self._tmp_dir) if os.path.splitext(filename)[-1] in ext]
		return [os.path.join(self._tmp_dir, filename) for filename in file_list]

	# Parse NUCmer delta file to get total alignment length and total sim_errors
	@staticmethod
	def parse_delta(filename):
		""" Reads a NUCmer output .delta file, extracting the aligned length and
			number of similarity errors for each aligned uniquely-matched region,
			and returns the cumulative total for each as a tuple.

			- filename is the path to the input .delta file
		"""
		aln_length, sim_errors = 0, 0
		for line in [l.strip().split() for l in open(filename, 'rU').readlines()]:
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
		""" Returns last exception as a string, or use in logging.
		"""
		exc_type, exc_value, exc_traceback = sys.exc_info()
		return ''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))

	# Parse NUCmer delta output to store alignment total length, sim_error,
	# and percentage identity, for each pairwise comparison
	def process_delta(self):
		""" Returns a tuple containing a list and four dictionaries. The list
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
		"""
		delta_files = self.get_input_files('.delta')
		self._logger.info("Delta files:\n\t%s" % '\n\t'.join(delta_files))
		self._logger.info("Processing .delta files")
		# We store pairwise comparison lengths in dictionaries, keyed by organism
		# ID pair tuples:
		# perc_aln is useful, as it is a matrix of the minimum percentage of an
		# organism's genome involved in a pairwise alignment
		lengths, sim_errors, perc_ids, perc_aln = {}, {}, {}, {}
		for delta_filename in delta_files:
			self._logger.info("Processing %s" % delta_filename)
			# TODO: not from filename!!
			qname, sname = os.path.splitext(os.path.split(delta_filename)[-1])[0].split('_vs_')
			self._logger.info("Query organism: %s; Subject organism: %s" % (qname, sname))
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
					self._logger.error("Please retry the NUCmer comparison for %s vs %s manually" % (qname, sname))
				else:
					self._logger.error(
						"The NUCmer comparison between %s and %s " % (qname, sname) +
						"has no usable output. The comparison may be " +
						"too distant for use. Consider using --maxmatch.")
				self._logger.error(self.last_exception())
			candidate_id = qname
			reference_id = sname
			if qname not in self._query_gid_to_location:
				candidate_id = sname
				reference_id = qname
			lengths[(candidate_id, reference_id)] = tot_length
			sim_errors[(candidate_id, reference_id)] = tot_sim_error
			perc_ids[(candidate_id, reference_id)] = perc_id
			perc_aln[(candidate_id, reference_id)] = 1. * tot_length / self._total_lengths[candidate_id]
			# perc_aln[sname] = 1.*tot_length/self.total_lengths[sname]
		return lengths, sim_errors, perc_ids, perc_aln

	# METHOD: ANIm
	# This method uses NUCmer to calculate pairwise alignments for the input
	# organisms, without chopping sequences into fragments. We follow the method
	# of Richter et al. (2009)
	def calculate_anim(self):
		""" Calculate ANI by the ANIm method, as described in Richter et al (2009)
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
		"""
		self._logger.info("Running ANIm method")
		self.run_cmd_lines()
		# lengths, sim_errors, perc_ids, perc_aln = self.process_delta()
		return self.process_delta()

	def calculate_best_anim(self):
		"""
		"""
		min_lengths, min_sim_errors, min_perc_ids, min_perc_aln, ncbi, best_ani_by_candidate_id = {}, {}, {}, {}, {}, {}
		lengths, sim_errors, perc_ids, perc_aln = self.calculate_anim()

		for candidate_id_reference_id in perc_ids:
			candidate_id = candidate_id_reference_id[0]
			reference_id = candidate_id_reference_id[1]
			# self._logger.info("{}: candidate_id: {}; reference_id: {}".format(candidate_id_reference_id, candidate_id, reference_id))
			# ani_ish = perc_ids[candidate_id_reference_id] * perc_aln[candidate_id_reference_id]
			if candidate_id not in best_ani_by_candidate_id or (perc_aln[candidate_id_reference_id] > self._minimum_alignment and best_ani_by_candidate_id[candidate_id] < perc_ids[candidate_id_reference_id]):
				best_ani_by_candidate_id[candidate_id] = perc_ids[candidate_id_reference_id]
				min_lengths[candidate_id] = lengths[candidate_id_reference_id]
				min_sim_errors[candidate_id] = sim_errors[candidate_id_reference_id]
				min_perc_ids[candidate_id] = perc_ids[candidate_id_reference_id]
				min_perc_aln[candidate_id] = perc_aln[candidate_id_reference_id]
				ncbi[candidate_id] = self._reference_gid_to_taxid[reference_id]  # reference_id.split('.')[0]
		return min_lengths, min_sim_errors, min_perc_ids, min_perc_aln, ncbi
