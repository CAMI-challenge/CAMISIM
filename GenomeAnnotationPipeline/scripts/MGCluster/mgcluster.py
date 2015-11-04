__author__ = 'hofmann'


import os
import glob
import time
import shutil
import tempfile
from scripts.Validator.validator import Validator


class MGCluster(Validator):
	"""
	Alignment and clustering of marker genes with references
	"""

	_label = "MGCluster"

	_cluster_method_choices = ['average', 'furthest', 'nearest']
	_file_name_map = "map.tsv"
	_silva_ref_files = ["mothur_ref_distances", "mothur_ref_names", "mothur_alignment_ref.fasta", _file_name_map]

	# mothur command: cluster.split
	_mothur_cmd_cluster_split = "; ".join([
		'unique.seqs(fasta={mg_fasta})',
		'align.seqs(candidate=current, template={ref_align}, align=gotoh, flip=t, processors={processors})',
		'remove.seqs(accnos={filename}.unique.flip.accnos, fasta=current, name=current)',
		'merge.files(input={filename}.names-{ref_names}, output={filename}.merged.names)',
		'merge.files(input={filename}.pick.names-{ref_names}, output={filename}.merged.names)',
		'set.current(name={filename}.merged.names, column={local_dist})',
		'dist.seqs(oldfasta={ref_align}, column=current, cutoff={cutoff}, processors={processors}, calc=onegap, countends=F)',
		'set.current(column={local_dist})',
		'cluster.split(cutoff={cutoff}, method={method}, precision={precision}, column={local_dist}, name={filename}.merged.names)'
	])

	# mothur command: cluster
	_mothur_cmd_cluster = "; ".join([
		"unique.seqs(fasta={mg_fasta})",
		"align.seqs(candidate=current, template={ref_align}, align=gotoh, flip=t, processors={processors})",
		"remove.seqs(accnos={filename}.unique.flip.accnos, fasta=current, name=current)",
		"merge.files(input={filename}.names-{ref_names}, output={filename}.merged.names)",
		"merge.files(input={filename}.pick.names-{ref_names}, output={filename}.merged.names)",
		"set.current(name={filename}.merged.names, column={local_dist})",
		"dist.seqs(oldfasta={ref_align}, column=current, cutoff={cutoff}, processors={processors}, calc=onegap, countends=F)",
		"set.current(column={local_dist})",
		"cluster(cutoff={cutoff}, method={method}, precision={precision}, column={local_dist}, name={filename}.merged.names)"
	])

	def __init__(
		self, mothur_executable, directory_silva_reference, max_processors=1, temp_directory=None,
		logfile=None, verbose=False, debug=False):
		"""
		Constructor

		@param mothur_executable: File path to mothur binary
		@type mothur_executable: str | unicode
		@param directory_silva_reference: Path to directory with SILVA reference database files
		@type directory_silva_reference: str | unicode
		@param max_processors: Maximum number of available processors
		@type max_processors: int | long
		@param temp_directory: Directory for temporary data
		@type temp_directory: str | unicode
		@param logfile: file handler or file path to a log file
		@type logfile: file | FileIO | StringIO | basestring
		@param verbose: Not verbose means that only warnings and errors will be past to stream
		@type verbose: bool
		@param debug: Display debug messages
		@type debug: bool
		"""
		assert self.validate_file(mothur_executable, executable=True)
		assert self.validate_dir(directory_silva_reference, file_names=self._silva_ref_files)
		assert self.validate_number(max_processors, minimum=1)
		assert self.validate_dir(temp_directory)
		super(MGCluster, self).__init__(logfile=logfile, verbose=verbose, debug=False)
		self._mothur_executable = mothur_executable
		self._tmp_dir = tempfile.mkdtemp(dir=temp_directory)
		self._max_processors = max_processors
		self._debug = debug
		ref_silva_distances = self._get_symbolic_link_path(os.path.join(directory_silva_reference, "mothur_ref_distances"))
		ref_silva_names = self._get_symbolic_link_path(os.path.join(directory_silva_reference, "mothur_ref_names"))  # unique
		ref_silva_alignment = self._get_symbolic_link_path(os.path.join(directory_silva_reference, "mothur_alignment_ref.fasta"))
		self._ref_silva_distances = ref_silva_distances
		self._ref_silva_names = ref_silva_names
		self._ref_silva_alignment = ref_silva_alignment
		# local_distance = os.path.join(self._working_dir, "ref.align.dist")
		self._local_distance = "ref.align.dist"

	def __exit__(self, type, value, traceback):
		super(MGCluster, self).__exit__(type, value, traceback)
		if not self._debug:
			shutil.rmtree(self._tmp_dir)

	def cluster(self, marker_gene_fasta, output_cluster_file, distance_cutoff, precision=1000, method="average"):
		"""
		CLuster Markergenes

		@param marker_gene_fasta: Fasta formatted file with marker genes
		@type marker_gene_fasta: str | unicode
		@param output_cluster_file: Output of mg clustering
		@type output_cluster_file: str | unicode
		@param distance_cutoff: Exclude irrelevant higher distances before clustering
		@type distance_cutoff: int | long
		@param precision: Cluster are made in steps: 10: 0.1, 100: 0.01, 1000: 0.001
		@type precision: int | long
		@param method: Cluster algorithm 'average', 'furthest', 'nearest'
		@type method: str | unicode

		@rtype: None
		"""
		assert self.validate_file(marker_gene_fasta)
		assert self.validate_dir(output_cluster_file, only_parent=True)
		assert self.validate_number(distance_cutoff, minimum=0, maximum=1)
		assert self.validate_number(precision, minimum=0)
		assert method in self._cluster_method_choices
		self._logger.info("Starting clustering process")
		start = time.time()
		old_dir = os.getcwd()
		# local paths required or mothur messes up the dist.seqs command, do not use absolut paths!!!
		os.chdir(self._tmp_dir)
		local_marker_gene_fasta = self._get_symbolic_link_path(marker_gene_fasta)
		shutil.copy2(self._ref_silva_distances, self._local_distance)

		mothur_cmd = self._get_mothur_cmd(local_marker_gene_fasta, distance_cutoff, precision, method=method)
		cmd = "{mothur_executable} '#{mothur_cmd}'".format(
			mothur_executable=self._mothur_executable,
			mothur_cmd=mothur_cmd)
		os.system(cmd)
		os.chdir(old_dir)

		project_folder = os.path.dirname(output_cluster_file)
		find_mask_list = os.path.join(self._tmp_dir, "*.list")
		list_of_files = glob.glob(find_mask_list)
		if len(list_of_files) == 0:
			msg = "Clustering with mothur failed #1"
			self._logger.error(msg)
			raise RuntimeError(msg)
		elif len(list_of_files) == 1:
			local_distance = os.path.join(self._tmp_dir, "ref.align.dist")
			if os.path.exists(local_distance):
				if self._debug:
					shutil.copy2(local_distance, os.path.join(project_folder, "mothur_distances.tsv"))
				shutil.copy2(list_of_files[0], output_cluster_file)
				self._logger.info("Clustering success")
			else:
				msg = "Clustering with mothur failed #2"
				self._logger.error(msg)
				raise RuntimeError(msg)
		else:
			msg = "Clustering with odd result, several files found!"
			self._logger.error(msg)
			raise RuntimeError(msg)
		end = time.time()

		# move logfiles
		find_mask_list = os.path.join(self._tmp_dir, "*.logfile")
		list_of_log_files = glob.glob(find_mask_list)
		for log_file in list_of_log_files:
			log_file_name = os.path.basename(log_file)
			shutil.copy2(log_file, os.path.join(project_folder, log_file_name))
		self._logger.info("Done ({}s)".format(round(end - start), 1))

	def _get_symbolic_link_path(self, original_file_path):
		"""
		Get path to local symbolic link since mothur might act odd else.

		@param original_file_path:
		@type original_file_path: str | unicode

		@return: Local path
		@rtype: str | unicode
		"""
		assert isinstance(original_file_path, basestring)
		basename = os.path.basename(original_file_path)
		new_path = os.path.join(self._tmp_dir, basename)
		os.symlink(original_file_path, new_path)
		# return new_path
		return basename

	def _get_mothur_cmd(self, marker_gene_fasta, cutoff, precision, method="average"):
		"""
		Get command line to run mothur

		@param marker_gene_fasta: Fasta formatted file with marker genes
		@type marker_gene_fasta: str | unicode
		@param cutoff: Exclude irrelevant higher distances before clustering
		@type cutoff: int | long
		@param precision: Cluster are made in steps: 10: 0.1, 100: 0.01, 1000: 0.001
		@type precision: int | long
		@param method: Cluster algorithm 'average', 'furthest', 'nearest'
		@type method: str | unicode

		@return: Command line
		@rtype: str | unicode
		"""
		assert self.validate_file(marker_gene_fasta)
		assert self.validate_number(cutoff, minimum=0, maximum=1)
		assert self.validate_number(precision, minimum=0)
		assert method in self._cluster_method_choices
		# basename = os.path.basename(marker_gene_fasta)
		# filename, extension = os.path.splitext(basename)
		filename, extension = os.path.splitext(marker_gene_fasta)
		#
		# mothur_cmd = MGCluster._mothur_cmd_cluster
		mothur_cmd = MGCluster._mothur_cmd_cluster_split
		return mothur_cmd.format(
			wd=self._tmp_dir,
			debug=self._tmp_dir,
			# filename=os.path.join(self._working_dir, filename),
			filename=filename,
			mg_fasta=marker_gene_fasta,
			ref_align=self._ref_silva_alignment,
			ref_names=self._ref_silva_names,
			local_dist=self._local_distance,
			processors=self._max_processors,
			cutoff=cutoff,
			precision=precision,
			method=method)

	@staticmethod
	def get_file_name_of_map():
		return MGCluster._file_name_map
