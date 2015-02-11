__author__ = 'hofmann'


import os
import glob
import time
import shutil
import tempfile
from source.argumenthandler import ArgumentHandler
from source.logger import Logger


def main(options):
	assert isinstance(options, ArgumentHandler)
	mg_cluster = MGCluster(mothur_executable=options.binary_mothur,
						   directory_silva_reference=options.silva_reference_directory,
						   max_processors=options.processors,
						   temp_directory=options.temp_directory,
						   debug=options._debug_mode)

	return mg_cluster.cluster(marker_gene_fasta=os.path.join(options.project_directory, options.file_mg_16s),
						output_cluster_file=os.path.join(options.project_directory, options.file_cluster_mg_16s),
						distance_cutoff=options.distance_cutoff,
						precision=options.precision,
						method=options.cluster_method)


class MGCluster(object):

	_mothur_cmd_ref_dist_split = """unique.seqs(fasta={mg_fasta})
align.seqs(candidate=current, template={ref_align}, align=gotoh, flip=t, processors={processors})
remove.seqs(accnos={filename}.unique.flip.accnos, fasta=current, name=current)
merge.files(input={filename}.names-{ref_names}, output={filename}.merged.names)
merge.files(input={filename}.pick.names-{ref_names}, output={filename}.merged.names)
set.current(name={filename}.merged.names, column={local_dist})
dist.seqs(oldfasta={ref_align}, column=current, cutoff={cutoff}, processors={processors}, calc=onegap, countends=F)
set.current(column={local_dist})
cluster.split(cutoff={cutoff}, method={method}, precision={precision}, column={local_dist}, name={filename}.merged.names)"""

	_mothur_cmd_ref_dist = """unique.seqs(fasta={mg_fasta})
align.seqs(candidate=current, template={ref_align}, align=gotoh, flip=t, processors={processors})
remove.seqs(accnos={filename}.unique.flip.accnos, fasta=current, name=current)
merge.files(input={filename}.names-{ref_names}, output={filename}.merged.names)
merge.files(input={filename}.pick.names-{ref_names}, output={filename}.merged.names)
set.current(name={filename}.merged.names, column={local_dist})
dist.seqs(oldfasta={ref_align}, column=current, cutoff={cutoff}, processors={processors}, calc=onegap, countends=F)
set.current(column={local_dist})
cluster(cutoff={cutoff}, method={method}, precision={precision}, name={filename}.merged.names)"""

	def __init__(self, mothur_executable, directory_silva_reference, max_processors=1, temp_directory=None, debug=False, logger=None):
		self._logger = logger
		if not self._logger:
			self._logger = Logger("MGCluster")
		self._mothur_executable = mothur_executable
		self._temp_directory = temp_directory
		self._working_dir = tempfile.mkdtemp(dir=self._temp_directory)

		self._old_dir = os.getcwd()
		# required or mothur messes up the dist.seqs command
		# do not use absolut paths!!!
		os.chdir(self._working_dir)

		self._max_processors = max_processors
		self._debug = debug
		ref_silva_distances = self._get_symbolic_link_path(os.path.join(directory_silva_reference, "mothur_ref_distances"))
		ref_silva_names = self._get_symbolic_link_path(os.path.join(directory_silva_reference, "mothur_ref_names"))  # unique
		ref_silva_alignment = self._get_symbolic_link_path(os.path.join(directory_silva_reference, "mothur_alignment_ref.fasta"))
		self._ref_silva_distances = ref_silva_distances
		self._ref_silva_names = ref_silva_names
		self._ref_silva_alignment = ref_silva_alignment
		#local_distance = os.path.join(self._working_dir, "ref.align.dist")
		self._local_distance = "ref.align.dist"

	def cluster(self, marker_gene_fasta, output_cluster_file, distance_cutoff, precision=1000, method="average"):
		self._logger.info("[MGCluster] Starting clustering process")
		start = time.time()
		local_marker_gene_fasta = self._get_symbolic_link_path(marker_gene_fasta)
		shutil.copy2(self._ref_silva_distances, self._local_distance)

		mothur_cmd = self._get_mothur_cmd(local_marker_gene_fasta, distance_cutoff, precision, method=method)
		cmd = "echo \"{mothur_cmd}\" | {mothur_executable}".format(mothur_cmd=mothur_cmd,
																   mothur_executable=self._mothur_executable)
		os.system(cmd)
		os.chdir(self._old_dir)

		find_mask_list = os.path.join(self._working_dir, "*.list")
		list_of_files = glob.glob(find_mask_list)
		if len(list_of_files) == 0:
			self._logger.error("[MGCluster] Clustering failed")
			self._logger.warning("[MGCluster] Remove manually: '{}'".format(self._working_dir))
			result = False
		elif len(list_of_files) == 1:
			local_distance = os.path.join(self._working_dir, "ref.align.dist")
			if os.path.exists(local_distance):
				parent_folder = os.path.dirname(output_cluster_file)
				shutil.copy2(local_distance, os.path.join(parent_folder, "mothur_distances.tsv"))
				shutil.copy2(list_of_files[0], output_cluster_file)
				self._logger.info("[MGCluster] Clustering success")
			else:
				self._logger.error("[MGCluster] Clustering failed!")
			if not self._debug:
				shutil.rmtree(self._working_dir)
			else:
				self._logger.warning("[MGCluster] Remove manually: '{}'".format(self._working_dir))
			result = True
		else:
			self._logger.error("[MGCluster] Clustering with odd result, several files found!")
			self._logger.warning("[MGCluster] Remove manually: '{}'".format(self._working_dir))
			result = False
		end = time.time()
		self._logger.info("[MGCluster] Done ({}s)".format(round(end - start), 1))
		return result

	def _get_symbolic_link_path(self, original_file_path):
		basename = os.path.basename(original_file_path)
		new_path = os.path.join(self._working_dir, basename)
		os.symlink(original_file_path, new_path)
		#return new_path
		return basename

	def _get_mothur_cmd(self, marker_gene_fasta, cutoff, precision, method="average"):
		#basename = os.path.basename(marker_gene_fasta)
		#filename, extension = os.path.splitext(basename)
		filename, extension = os.path.splitext(marker_gene_fasta)
		#
		#mothur_cmd = MGCluster._mothur_cmd_ref_dist
		mothur_cmd = MGCluster._mothur_cmd_ref_dist_split
		return mothur_cmd.format(wd=self._working_dir,
								 debug=self._working_dir,
								 #filename=os.path.join(self._working_dir, filename),
								 filename=filename,
								 mg_fasta=marker_gene_fasta,
								 ref_align=self._ref_silva_alignment,
								 ref_names=self._ref_silva_names,
								 local_dist=self._local_distance,
								 processors=self._max_processors,
								 cutoff=cutoff,
								 precision=precision,
								 method=method)
